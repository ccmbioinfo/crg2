if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/{contig}.regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


def get_decoy_removed_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("decoy_rm/{family}-{sample}.no_decoy_reads.bam",
                  sample=wildcards.sample,
                  family=wildcards.family)

rule drag_str:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        dragstr_table=config["ref"]["dragstr_table"]
    output:
        dragstr_parameters="called/{family}_{sample}.drag_str_parameters.txt"
    log:
        "logs/gatk/CalibrateDragstrModel/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "CalibrateDragstrModel")


if config["run"]["gatk"] == "gatk":
    rule gatk:
        input:
            bam=get_sample_bams,
            ref=config["ref"]["genome"],
            known=config["ref"]["known-variants"],
            regions="called/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
        output:
            gvcf=temp("called/{family}_{sample}.{contig}.g.vcf.gz")
        log:
            "logs/gatk/haplotypecaller/{family}_{sample}.{contig}.log"
        params:
            extra=get_call_variants_params,
            java_opts=config["params"]["gatk"]["java_opts"],
        group: "gatkcall"
        resources:
            mem=lambda wildcards, input: len(input.bam) * 15
        wrapper:
            get_wrapper_path("gatk", "haplotypecaller")
elif config["run"]["gatk"] == "dragen":
    rule dragen:
        input:
            bam=get_sample_bams,
            ref=config["ref"]["genome"],
            known=config["ref"]["known-variants"],
            regions="called/{contig}.regions.bed" if config["processing"].get("restrict-regions") else [],
            dragstr_parameters="called/{family}_{sample}.drag_str_parameters.txt"
        output:
            gvcf=temp("called/{family}_{sample}.{contig}.g.vcf.gz")
        log:
            "logs/gatk/dragenhaplotypecaller/{family}_{sample}.{contig}.log"
        params:
            extra=get_call_variants_params,
            java_opts=config["params"]["dragen"]["java_opts"],
        group: "gatkcall"
        resources:
            mem=lambda wildcards, input: len(input.bam) * 15
        wrapper:
            get_wrapper_path("gatk", "dragenhaplotypecaller") 


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/{{family}}_{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf=temp("called/{family}.{contig}.g.vcf.gz")
    params:
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/combinegvcfs.{family}.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "combinegvcfs")


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/{family}.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/{family}.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/genotypegvcfs.{family}.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "genotypegvcfs")


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("genotyped/{family}.{contig}.vcf.gz", contig=get_canon_contigs(), family=project),
	## use this to remove repetitive contigs for dag generation
	#vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig="GRCh37"), 
    output:
        vcf=protected("genotyped/{family}.vcf.gz")
    log:
        "logs/picard/{family}.merge-genotyped.log"
    wrapper:
        get_wrapper_path("picard", "mergevcfs")
