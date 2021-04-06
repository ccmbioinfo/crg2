rule vep:
    input:
        "annotated/{project}-ensemble-decomposed.vcf.gz",
    output:
        #temp("annotated/vep/{project}-ensemble-decomposed.vep.vcf"),
        "annotated/vep/{project}-ensemble-decomposed.vep.vcf",
    log:
        "logs/vep/vep.{project}.log"
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
        "annotated/vep/{project}-ensemble-decomposed.vep.vcf"
    output:
        "annotated/vcfanno/{project}-ensemble-decomposed.vep.vcfanno.vcf"
    log:
        "logs/vcfanno/vcfanno.{project}.log"
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
        "annotated/vcfanno/{project}-ensemble-decomposed.vep.vcfanno.vcf"
    output:
         db = "annotated/{project}-gemini.db",
    log:
        "logs/vcf2db/vcf2db.{project}.log"
    params:
        ped = config["run"]["ped"],
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")

