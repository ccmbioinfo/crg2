def get_cre_bams(ext="bam"):
    if gatk == "gatk3":
        return expand("recal/gatk3/{family}_{sample}.{ext}", family=project, sample=samples.index, ext=ext)
    return expand("recal/{family}_{sample}.{ext}", family=project, sample=samples.index, ext=ext)

rule gathervcf:
    input:
        vcfs = lambda w: expand("called/gatk3/{family}-{contig}.vcf", family=project, contig=get_canon_contigs()),
    output:
        gvcf=protected("genotyped/{family}-gatk3_haplotype.vcf")
    wrapper:
        get_wrapper_path("picard","gathervcfs")


#duplicating gatk4 rules from crg2/calling.smk for cre file namings 
#sub-workflows does not seem to work smoothly
rule call_variants:
    input:
        bam=get_sample_bams,
        #bam=get_cre_bams(),
        ref=config["ref"]["genome"],
        known=config["ref"]["known_variants"],
        regions="mapped/bed/{family}-sort-callable-{contig}.bed",
        #regions="called/gatk/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf="called/gatk/{family}_{sample}.{contig}.g.vcf.gz"
    log:
        "logs/gatk/haplotypecaller/{family}_{sample}.{contig}.log"
    params:
        extra=get_call_variants_params,
        java_opts=config["params"]["gatk"]["java_opts"],
    group: "gatkcall"
    resources: 
        mem=lambda wildcards, input: len(input.bam) * 20
    wrapper:
        get_wrapper_path("gatk", "haplotypecaller")


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/gatk/{{family}}_{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf=temp("called/gatk/{family}.{contig}.g.vcf.gz")
    params:
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/{family}.combinegvcfs.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "combinegvcfs")


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/gatk/{family}.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/gatk/{family}.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/{family}.genotypegvcfs.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "genotypegvcfs")


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("genotyped/gatk/{{family}}.{contig}.vcf.gz", contig=get_canon_contigs()),
	## use this to remove repetitive contigs for dag generation
	#vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig="GRCh37"), 
    output:
        vcf="genotyped/gatk/{family}.vcf.gz"
    log:
        "logs/picard/{family}-merge-genotyped.log"
    wrapper:
        get_wrapper_path("picard", "mergevcfs")

rule gatk4:
    input: 
        "genotyped/gatk/{family}.vcf.gz".format(family=project)
    output:
        gvcf=protected("genotyped/{family}-gatk_haplotype.vcf")
    log:
        "logs/gatk/{family}.log"
    shell:
        '''
        gunzip -c -d {input} > {output}
        '''
 
rule bgzip:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    conda:
        "../../envs/common.yaml"

    shell:
        '''
        bgzip -c {input} > {output}
        '''

rule tabix:
    input: 
        "{prefix}.vcf.gz"
    output: 
        "{prefix}.vcf.gz.tbi"
    log: 
        "logs/{prefix}.log"
    conda:
           "../../envs/common.yaml"
    wrapper:
        get_wrapper_path("tabix")
