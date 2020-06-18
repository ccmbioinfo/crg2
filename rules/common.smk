import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.7.1")

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["run"]["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str)
                                      for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

project = config["run"]["project"]
flank = config["run"]["flank"]

##### Wildcard constraints #####
wildcard_constraints:
    vartype = "snvs|indels",
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"])


##### Helper functions #####

def get_fai():
    return config["ref"]["genome"] + ".fai"


# contigs in reference genome
def get_contigs():
    return pd.read_table(get_fai(),
                         header=None, usecols=[0], squeeze=True, dtype=str)


def is_autosomal(chrom):
    # from bcbio-nextgen/bcbio/heterogeneity/chromhacks.py
    """Keep chromosomes that are a digit 1-22, or chr prefixed digit chr1-chr22
    """
    try:
        int(chrom)
        return True
    except ValueError:
        try:
            int(str(chrom.lower().replace("chr", "").replace("_", "").replace("-", "")))
            return True
        except ValueError:
            return False


def is_sex(chrom):
    return chrom in ["X", "chrX", "Y", "chrY"]


def is_mitochondrial(chrom):
    return chrom in ["MT", "chrM", "chrMT"]


def is_autosomal_or_sex(chrom):
    return is_autosomal(chrom) or is_sex(chrom)


def is_nonalt(chrom):
    """Check that a chromosome is on 1-22, X, Y, MT.
    """
    return is_autosomal_or_sex(chrom) or is_mitochondrial(chrom)


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    sample = wildcards.sample
    unit = wildcards.unit
    fastqs = units.loc[(sample, unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return [fastqs.fq1, fastqs.fq2]
    elif len(fastqs) == 1:
        return [fastqs.fq1]
    else:
        return ["fastq/{sample}-{unit}_1.fq".format(sample=sample, unit=unit), "fastq/{sample}-{unit}_2.fq".format(sample=sample, unit=unit)]


def get_bam(wildcards):
    """Get previously aligned bam files of given sample-unit."""
    bam_file = units.loc[(wildcards.sample, wildcards.unit), ["bam"]].dropna()
    return bam_file.bam


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


def get_wrapper_path(*dirs):
    return "file:%s" % os.path.join(workflow.basedir, "wrappers", *dirs)
