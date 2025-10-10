def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="genotyped/gatk/{family}.vcf.gz"
    output:
        vcf=temp("filtered/{family}.{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{family}.{vartype}.log"
    wrapper:
        get_wrapper_path("gatk", "selectvariants")


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="filtered/{family}.{vartype}.vcf.gz"
    output:
        vcf=temp("filtered/{family}.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        "logs/gatk/variantfiltration/{family}.{vartype}.log"
    wrapper:
        get_wrapper_path("gatk", "variantfiltration")


rule recalibrate_calls:
    input:
        vcf="filtered/{family}.{vartype}.vcf.gz"
    output:
        vcf=temp("filtered{family}.{vartype}.recalibrated.vcf.gz")
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        "logs/gatk/variantrecalibrator/{family}.{vartype}.log"
    wrapper:
        get_wrapper_path("gatk", "variantrecalibrator")


rule merge_calls:
    input:
        vcfs=expand("filtered/{family}.{vartype}.{filtertype}.vcf.gz",
                   vartype=["snvs", "indels"],
                   family=project,
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf="filtered/{family}.gatk_haplotype.vcf.gz"
    log:
        "logs/picard/{family}.merge-filtered.log"
    wrapper:
        get_wrapper_path("picard", "mergevcfs")

