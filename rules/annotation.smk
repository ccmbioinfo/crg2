rule vt:
    input:
        "filtered/all.vcf.gz", # (vcf, bcf, or vcf.gz)
    output:
        temp("filtered/all.uniq.normalized.decomposed.vcf"),
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/vt.uniq.normalized.decomposed.log"
    wrapper:
        get_wrapper_path("vt")

rule vep:
    input:
        "filtered/all.uniq.normalized.decomposed.vcf",
    output:
        temp("annotated/vep/all.vep.vcf"),
    log:
        "logs/vep/vep.log"
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
        "annotated/vep/all.vep.vcf",
    output:
        "annotated/vcfanno/all.vcfanno.vcf",
    log:
        "logs/vcfanno/vcfanno.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
        lua_script=config["annotation"]["vcfanno"]["lua_script"],
       	conf=config["annotation"]["vcfanno"]["conf"],
        base_path=config["annotation"]["vcfanno"]["base_path"],
    wrapper:
        get_wrapper_path("vcfanno")

rule vcf2db:
    input:
        "annotated/vcfanno/all.vcfanno.vcf",
    output:
         db="annotated/gemini.db",   # annotated calls (vcf, bcf, or vcf.gz)
    log:
        "logs/vcf2db/vcf2db.log"
    params:
        ped=config["annotation"]["vcf2db"]["ped"],
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")
