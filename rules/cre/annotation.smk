rule vep:
    input:
        "annotated/coding/{family}.ensemble.decomposed.vcf.gz",
    output:
        temp("annotated/coding/vep/{family}.coding.vep.vcf")
    log:
        "logs/vep/{family}.vep.coding.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        ref=config["ref"]["genome"],
    wrapper:
        get_wrapper_path("vep")

rule vcfanno:
    input:
        "annotated/coding/vep/{family}.coding.vep.vcf"
    output:
        "annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf"
    log:
        "logs/vcfanno/{family}.vcfanno.coding.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
        lua_script = config["annotation"]["cre.vcfanno"]["lua_script"],
       	conf = config["annotation"]["cre.vcfanno"]["conf"],
        base_path = config["annotation"]["cre.vcfanno"]["base_path"],
    wrapper:
        get_wrapper_path("vcfanno")



rule vcf2db:
    input:
        "annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf".format(family=project)
    output:
         db = "annotated/coding/{family}-gemini.db",
    log:
        "logs/vcf2db/vcf2db.{family}.log"
    params:
        ped = config["run"]["ped"],
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")


