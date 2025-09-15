rule vep_mosaic:
    input:
        "filtered/{family}-{p}.uniq.normalized.decomposed.vcf.gz"
    output:
        temp("annotated/{p}/vep/{family}.coding.vep.vcf")
    log:
        "logs/vep/{p}/{family}.vep.coding.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        ref=config["ref"]["genome"],
    wrapper:
        get_wrapper_path("vep")


