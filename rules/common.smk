
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
from datetime import date

min_version("5.7.1")

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["run"]["units"], dtype=str).set_index(["sample"], drop=False)

validate(units, schema="../schemas/units.schema.yaml")

project = config["run"]["project"]
flank = config["run"]["flank"]
gatk = config["run"]["gatk"] 

##### Wildcard constraints #####
wildcard_constraints:
    vartype = "snvs|indels",
    sample = "|".join(samples.index),
    family = project


##### Helper functions #####

def get_fai():
    return config["ref"]["genome"] + ".fai"


# contigs in reference genome
def get_contigs():
    return pd.read_table(get_fai(),
                         header=None, usecols=[0], squeeze=True, dtype=str)

def get_canon_contigs():
    return pd.read_table(config["ref"]["canon_bed"],
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


def get_bam(wildcards):
    """Get previously aligned bam files of given family-sample."""
    bam_file = units.loc[(wildcards.sample), ["bam"]].dropna()
    return bam_file.bam


def is_single_end(sample, unit):
    """Return True if family-sample is single end."""
    return pd.isnull(units.loc[(sample), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{family}_{sample}\tSM:{family}_{sample}\tPL:{platform}'".format(
        family=wildcards.family,
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample), "platform"])

def get_sample_order(wildcards):
    """bcftools view -s <sample_order> as freebaayes randomly picks order"""
    return [ "{}_{}".format(project,s) for s in samples.index ]

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given family-sample."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{family}_{sample}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{family}_{sample}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{family}_{sample}.bam",
                  sample=wildcards.sample,
                  family=wildcards.family)


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {} ".format(wildcards.contig)) +
            config["params"][gatk]["HaplotypeCaller"])


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{family}_{sample}.sorted.bam"
    if config["processing"]["mark-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{family}_{sample}.bam"
    if bai or  config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
    return f
    

def get_recal_input_gatk3(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{family}_{sample}.sorted.bam"
    if config["processing"]["mark-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{family}_{sample}.bam"
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

def get_annotated_sv_vcf():
    """Get the annotated MetaSV vcf of given sample."""
    return ["sv/metasv/{family}_{sample}/variants.snpeff.svscore.vcf".format(sample=sample, family=project) for sample in samples.index]

def get_wrapper_path(*dirs):
    return "file:%s" % os.path.join(workflow.basedir, "wrappers", *dirs)

def get_eh_json(wildcards):
    """Get the EH JSON of all samples."""
    return ["str/EH/{}_{}.json".format(wildcards.family, sample) for sample in samples.index]

def parse_ped_id(individual_id, family):
    if individual_id != "0":
        parsed_id = family + "_" + individual_id.replace(family, "")
    else:
        parsed_id = "0"

    return parsed_id

def format_pedigree(wildcards):
    family = wildcards.family
    ped = config["run"]["ped"]
    if ped == "":
        return None
    else:
        ped = pd.read_csv(
            ped,
            sep=" ",
            header=None,
            names=["fam_id", "individual_id", "pat_id", "mat_id", "sex", "phenotype"],
        )

        ped["fam_id"] = family
        for col in ["individual_id", "pat_id", "mat_id"]:
            ped[col] = [parse_ped_id(individual_id, family) for individual_id in ped[col].values]

        ped.to_csv(f"{family}.ped", sep=" ", index=False, header=False)

        return f"{family}.ped"

