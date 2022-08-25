import pandas as pd
from datetime import datetime

configfile: "config.yaml"

samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)
project = config["run"]["project"]

rule mity_call:
    input:
        bam=[expand("recal/{family}_{sample}.bam".format(family=config["run"]["project"], sample=s)) for s in samples.index]
    output:
        protected("mitochondrial_variants/{family}.pre-normalised.mity.vcf.gz")
    params:
        outdir="mitochondrial_variants/",
        prefix="{family}.pre-normalised",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_call/{family}.mity_call.log"
    wrapper:    
        get_wrapper_path("mity/call")

rule mity_normalise:
    input:
        "mitochondrial_variants/{family}.pre-normalised.mity.vcf.gz"
    output:
        protected("mitochondrial_variants/{family}.normalised.mity.vcf.gz")
    params:
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_normalise/{family}.mity_normalise.log"
    wrapper:    
        get_wrapper_path("mity/normalise")

rule mito_vcfanno:
    input:
        "mitochondrial_variants/{family}.normalised.mity.vcf.gz"
    output:
        "mitochondrial_variants/vcfanno/{family}.normalised.mity.vcfanno.vcf"
    log:
        "logs/mity/vcfanno/{family}.mity.vcfanno.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
       	conf = config["annotation"]["mt.vcfanno"]["conf"],
        base_path = config["annotation"]["mt.vcfanno"]["base_path"]
    wrapper:
        get_wrapper_path("vcfanno")

rule mity_report:
    input:
        "mitochondrial_variants/vcfanno/{family}.normalised.mity.vcfanno.vcf"
    output:
        "report/mitochondrial/{family}.mitochondrial.report.csv"
    params:
        outdir="report/mitochondrial/",
        prefix="{family}_mito",
        tool=config["tools"]["mity"]
    log:
        "logs/mity/mity_report/{family}.mity_report.log"
    wrapper:
        get_wrapper_path("mity/report")