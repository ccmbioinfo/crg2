# Adapted from https://github.com/GooglingTheCancerGenome/sv-callers/blob/master/snakemake/rules/manta.smk
# Need to adapt this to handle multiple samples, i.e. outdir willl be sv/manta/<sample>
# Also don't want bam hard-coded
rule manta:  # single-sample analysis
    input:
        bam = "recal/NA12878-1.bam",
        fasta = config["ref"]["genome"]
    params:
        # excl_opt = '-x "%s"' % get_bed() if exclude_regions() else ""
        outdir = "sv/manta/"
    threads:
        4
    output:
        "sv/manta/results/variants/diploidSV.vcf.gz"
    wrapper:
        get_wrapper_path("manta")
