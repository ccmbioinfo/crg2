rule gatk_call_mosaic:
    input:
        map = "recal/{family}_{sample}.bam",
        fasta=config["ref"]["genome"],
        mgp_germline=config["ref"]["known_variants"]
    output:
        vcf="called/gatk_mutect/{family}_{sample}.vcf",
        stats="called/gatk_mutect/{family}_{sample}.vcf.stats"
    threads: 4
    log:
        "logs/gatk/mutect/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "mutect")


rule pileup_summaries:
    input:
        cram="recal/{family}_{sample}.cram",
        cram_index="recal/{family}_{sample}.cram.crai",
        common_variants=config["ref"]["known_variants_common"],
        ref=config["ref"]["genome"]
    output:
        "called/gatk_mutect/{family}_{sample}.pileups.table"
    log:
        "logs/gatk/getpileupsummaries/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "getpileupsummaries")


rule calculate_contamination:
    input:
        pileups="called/gatk_mutect/{family}_{sample}.pileups.table"
    output:
        segments="called/gatk_mutect/{family}_{sample}.segments.table",
        table="called/gatk_mutect/{family}_{sample}.contamination.table"
    log:
        "logs/gatk/calculatecontamination/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "calculatecontamination")
    

rule filter_mutect_call:
    input: 
        vcf="called/gatk_mutect/{family}_{sample}.vcf",
        fasta=config["ref"]["genome"],
        contamination="called/gatk_mutect/{family}_{sample}.contamination.table",
        segments="called/gatk_mutect/{family}_{sample}.segments.table"
    output:
        vcf="genotyped/{family}_{sample}-gatk_somatic.vcf"
    log:
        "logs/gatk/filtermutectcalls/{family}_{sample}_somatic.log"
    params:
        extra=config["params"]["gatk"]["FilterMutectCalls"]
        #java_opts=config["params"]["gatk"]["java_opts"]
    wrapper:
        get_wrapper_path("gatk", "filtermutectcalls")


rule pass_somatic:
    input:
        "genotyped/{family}_{sample}-gatk_somatic.vcf.gz", "genotyped/{family}_{sample}-gatk_somatic.vcf.gz.tbi"
    output:
        "genotyped/{family}_{sample}-gatk_somatic.pass_mutect2.vcf.gz"
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
        vcf_unsort=temp("genotyped/{family}-gatk_somatic.unsorted.vcf"),
        vcf="genotyped/{family}-gatk_somatic.vcf"
    log: 
        "logs/bcftools/merge/{family}_gatk_somatic.log"
    params:
        extra=config["params"]["bcftools"]["merge"]
    wrapper:
        get_wrapper_path("bcftools", "merge")
