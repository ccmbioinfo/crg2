rule run_mity:
    input:
        bam=[expand("recal/{family}_{sample}.bam".format(family=config["run"]["project"], sample=s)) for s in samples.index]
    output:
        "mitochondrial_variants/{family}.normalise.mity.annotated.vcf.gz",
        "mitochondrial_variants/{family}.annotated_variants.xlsx"
    params:
        outdir="mitochondrial_variants",
        prefix="{family}",
        tool=config["tools"]["mity"],
        vcfanno_config=config["annotation"]["mt.vcfanno"]["conf"]
    log:
        "logs/mity/{family}.mity_run_all.log"
    wrapper:
        get_wrapper_path("mity")

rule generate_mt_report:
    input:
        vcf="mitochondrial_variants/{family}.normalise.mity.annotated.vcf.gz",
        report="mitochondrial_variants/{family}.annotated_variants.xlsx"
    output:
        "report/mitochondrial/{family}.mitochondrial.report.csv"        
    log:
        "logs/report/mitochondrial/{family}.mitochondrial.report.log"
    conda:
        "../envs/mt_report.yaml"
    script:
        "../scripts/mt_report.py"