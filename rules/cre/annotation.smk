rule vep:
    input:
        "annotated/coding/{family}-ensemble-decomposed.vcf.gz",
    output:
        #temp("annotated/vep/{family}-ensemble-decomposed.vep.vcf"),
        "annotated/coding/vep/{family}-ensemble-decomposed.vep.vcf",
    log:
        "logs/vep/vep.{family}.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        maxentscan=config["annotation"]["vep"]["maxentscan"],
        human_ancestor_fasta=config["annotation"]["vep"]["human_ancestor_fasta"],
        ref=config["ref"]["genome"],
    wrapper:
        get_wrapper_path("vep")

rule vcfanno:
    input:
        "annotated/coding/vep/{family}-ensemble-decomposed.vep.vcf"
    output:
        "annotated/coding/vcfanno/{family}-ensemble-decomposed.vep.vcfanno.vcf"
    log:
        "logs/vcfanno/vcfanno.{family}.log"
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
        "annotated/coding/vcfanno/{family}-ensemble-decomposed.vep.vcfanno.vcf".format(family=project)
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


