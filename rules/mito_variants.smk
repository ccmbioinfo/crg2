
rule mity_call:
    input:
        bam=[expand("recal/{family}_{sample}.bam".format(family=config["run"]["project"], sample=s)) for s in samples.index]
    output:
        protected("mitochondrial_variants/{family}.mity.call.vcf.gz")
    params:
        outdir="mitochondrial_variants/",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_call/{family}.mity_call.log"
    wrapper:    
        get_wrapper_path("mity/call")

rule mity_normalise:
    input:
        "mitochondrial_variants/{family}.mity.call.vcf.gz"
    output:
        protected("mitochondrial_variants/{family}.mity.normalise.decompose.vcf.gz")
    params:
        outdir="mitochondrial_variants/",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_normalise/{family}.mity_normalise.log"
    wrapper:    
        get_wrapper_path("mity/normalise")

rule mity_report:
    input:
        "mitochondrial_variants/{family}.mity.normalise.decompose.vcf.gz"
    output:
        "mitochondrial_variants/{family}.mity.annotated.vcf",
        "mitochondrial_variants/{family}.annotated_variants.xlsx"
    params:
        outdir="mitochondrial_variants/",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_report/{family}.mity_report.log"
    wrapper:
        get_wrapper_path("mity/report")