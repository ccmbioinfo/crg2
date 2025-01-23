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