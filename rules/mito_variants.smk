import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)
project = config["run"]["project"]

rule mity_call:
    input:
        bam=[expand("recal/{family}_{sample}.bam".format(family=config["run"]["project"], sample=s)) for s in samples.index]
    output:
        protected("mitochondrial_variants/mito_{family}.mity.vcf.gz")
    params:
        outdir="mitochondrial_variants/",
        prefix="mito_{family}",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_call/{family}.mity_call.log"
    wrapper:    
        get_wrapper_path("mity/call")


rule mito_vcfanno:
    input:
        "mitochondrial_variants/mito_{family}.mity.vcf.gz"
    output:
        "mitochondrial_variants/vcfanno/{family}_mito.vcfanno.vcf"
    log:
        "logs/mity/vcfanno/{family}.mito.vcfanno.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
       	conf ="~/crg2/vcfanno/mt.vcfanno.conf",
        base_path = "/hpf/largeprojects/ccmbio/ajain/mity/vcfanno"
    wrapper:
        get_wrapper_path("vcfanno")

rule mity_report:
    input:
        "mitochondrial_variants/vcfanno/{family}_mito.vcfanno.vcf"
    output:
        "report/mitochondrial/{family}_mity_vcfanno_final_report.csv"
    params:
        outdir="report/mitochondrial/",
        prefix="{family}_mito",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_report/{family}.mity_report.log"
    wrapper:
        get_wrapper_path("mity/report")