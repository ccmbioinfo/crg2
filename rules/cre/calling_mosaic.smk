rule gatk_call_mosaic:
    input:
        map=get_sample_bams
        fasta=config["ref"]["genome"],
    output:
        vcf="called/gatk/mutect/{family}_{sample}.vcf.gz"
        #vcf=temp("called/gatk/{family}_{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/mutect/{family}_{sample}.log"
    params:
        extra=config["params"]["gatk"]["Mutect"]
        java_opts=config["params"]["gatk"]["java_opts"]

    group: "gatkcall"
    resources: 
        mem=lambda wildcards, input: len(input.bam) * 15
    wrapper:
        get_wrapper_path("gatk", "mutect")



