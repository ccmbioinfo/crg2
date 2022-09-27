

rule freebayes_mosaic:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        regions="mapped/bed/{family}-sort-callable-{{contig}}.bed".format(family=project)
    output:
        temp("called/freebayes_mosaic/{contig}.vcf") 
    log:
        "logs/freebayes/mosaic/{contig}.log"
    params:
        extra=config["params"]["freebayes"]["mosaic_call"], 
        chunksize=100000
    threads: 1
    resources:
        mem=lambda wildcards, threads: threads * 4 if threads > 1 else 20
    wrapper:
        get_wrapper_path("freebayes")
    

rule merge_freebayes_mosaic:
    input:
        expand("called/freebayes_mosaic/{contig}.vcf",contig = get_canon_contigs())
    output:
        protected("genotyped/{family}-freebayes_mosaic.vcf")
    log:
        "logs/freebayes/mosaic/{family}-merge.log"
    shell:
        '''
        bcftools concat {input} | bcftools sort > {output}  
        '''


rule gatk_call_mosaic:
    input:
        map = "recal/{family}_{sample}.bam",
        #map="genes/{family}_{sample}_slice.bam",
        fasta=config["ref"]["genome"],
        gnomad_germline=config["params"]["gatk"]["Mutect2"]["gnomad_germline"]
    output:
        vcf="called/gatk_mutect/{family}_{sample}.vcf",
        stats="called/gatk_mutect/{family}_{sample}.vcf.stats"
    log:
        "logs/gatk/mutect/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "mutect")


rule filter_mutect_call:
    input: 
        vcf="called/gatk_mutect/{family}_{sample}.vcf",
        fasta=config["ref"]["genome"]
    output:
        vcf="genotyped/{family}_{sample}-gatk_somatic.vcf"
    log:
        "logs/gatk/filtermutectcalls/{family}_{sample}_somatic.log"
    params:
        extra=config["params"]["gatk"]["FilterMutectCalls"]
        #java_opts=config["params"]["gatk"]["java_opts"]
    wrapper:
        get_wrapper_path("gatk", "filtermutectcalls")


# filtering: skip softfilter step to do rule pass. Changed file name to distinguish from  "uniq.normalized.decomposed.pass.vcf.gz"
rule pass_mosaic:
    input:
        "filtered/{family}_{sample}-gatk_somatic.uniq.normalized.decomposed.vcf.gz", "filtered/{family}_{sample}-gatk_somatic.uniq.normalized.decomposed.vcf.gz.tbi"
    output:
        "filtered/{family}_{sample}-gatk_somatic.uniq.normalized.decomposed.pass.mosaic.vcf.gz"
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: 
        samples = get_sample_order,
        filter = "-f 'PASS,.' "
    wrapper:
        get_wrapper_path("bcftools", "view")


rule merge_mutect_sample:
    input:
        vcf=get_gatk_somatic_vcf(),
        index=get_gatk_somatic_vcf(ext="vcf.gz.tbi")
    output:
        vcf_unsort=temp("filtered/{family}-gatk_somatic.uniq.normalized.decomposed.unsorted.pass.mosaic.vcf"),
        vcf="filtered/{family}-gatk_somatic.uniq.normalized.decomposed.pass.mosaic.vcf"
    log: 
        "logs/bcftools/merge/{family}_gatk_somatic.log"
    params:
        extra=config["params"]["bcftools"]["merge"]
    wrapper:
        get_wrapper_path("bcftools", "merge")
