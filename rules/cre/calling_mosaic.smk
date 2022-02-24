
rule gatk_call_mosaic:
    input:
        map=get_sample_bams,
        fasta=config["ref"]["genome"]
    output:
        vcf="called/gatk_mutect/{family}_{sample}.vcf",
        stats="called/gatk_mutect/{family}_{sample}.vcf.stats"
        #vcf=temp("called/gatk/{family}_{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/mutect/{family}_{sample}.log"
    params:
        extra=config["params"]["gatk"]["Mutect"],
        #java_opts=config["params"]["gatk"]["java_opts"]
    wrapper:
        get_wrapper_path("gatk", "mutect")


rule filter_mutect_call:
    input: 
        vcf="called/gatk_mutect/{family}_{sample}.vcf",
        fasta=config["ref"]["genome"]
    output:
        vcf="genotyped/gatk_mutect/{family}_{sample}_somatic.vcf"
    log:
        "logs/gatk/filtermutectcalls/{family}_{sample}_somatic.log"
    params:
        extra=config["params"]["gatk"]["FilterMutectCalls"]
        #java_opts=config["params"]["gatk"]["java_opts"]
    wrapper:
        get_wrapper_path("gatk", "filtermutectcalls")

rule merge_mutect_sample:
    input:
        vcf=get_gatk_somatic_vcf(),
        index=get_gatk_somatic_vcf(ext="vcf.gz.tbi")
    output:
        vcf_unsort="genotyped/{family}-gatk_somatic_unsorted.vcf",
        vcf="genotyped/{family}-gatk_somatic.vcf"
    log: 
        "logs/bcftools/merge/{family}_gatk_somatic.log"
    params:
        extra=config["params"]["bcftools"]["merge"]
    wrapper:
        get_wrapper_path("bcftools", "merge")

rule freebayes_mosaic:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        regions="mapped/bed/{family}-sort-callable-{{contig}}.bed".format(family=project)
        #regions=config["ref"]["canon_bed"]
    output:
        temp("called/freebayes_mosaic/{contig}.vcf")  # either .vcf or .bcf
    log:
        "logs/freebayes/mosaic/{contig}.log"
    params:
        extra=config["params"]["freebayes"]["mosaic_call"], #["mosaic_call"],         # optional parameters
        chunksize=100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 1
    resources:
        mem=lambda wildcards, threads: threads * 4 if threads > 1 else 20
    wrapper:
        get_wrapper_path("freebayes")
    
rule merge_freebayes_mosaic:
    input:
        expand("called/freebayes_mosaic/{contig}.vcf",contig = get_canon_contigs())
    output:
        protected("genotyped/{family}-freebayes-mosaic.vcf")
    log:
        "logs/freebayes/mosaic/{family}-merge.log"
    shell:
        '''
        bcftools concat {input} | bcftools sort > {output}  
        '''

# calls will be bgzipped and tabixed with rule bgzip and tabix from cre/calling.smk
