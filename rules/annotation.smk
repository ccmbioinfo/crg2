def get_filt_vcf(wildcards):
    if wildcards.p == "coding":
        return "filtered/all.vcf.gz"
    else:
        return "filtered/{}/all.{}.vcf.gz".format(wildcards.p,wildcards.p)


rule vt:
    input: get_filt_vcf # (vcf, bcf, or vcf.gz)
    output:
        temp("filtered/{p}/all.{p}.uniq.normalized.decomposed.vcf"),
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/vt.{p}.uniq.normalized.decomposed.log"
    wrapper:
        get_wrapper_path("vt")

rule pass:
    input:
       	"{prefix}.{ext}"
    output:
        temp("{prefix}.pass.{ext,(vcf|vcf\.gz)}")
    threads: 1
    resources:
        mem_mb = 4000
    params: "-f PASS"
    wrapper:
        get_wrapper_path("bcftools", "view")

rule vep:
    input:
        "filtered/{p}/all.{p}.uniq.normalized.decomposed.pass.vcf",
    output:
        temp("annotated/{p}/vep/all.{p}.vep.vcf"),
    log:
        "logs/vep/vep.{p}.log"
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
        "annotated/{p}/vep/all.{p}.vep.vcf",
    output:
        "annotated/{p}/vcfanno/all.{p}.vep.vcfanno.vcf",
    log:
        "logs/vcfanno/vcfanno.{p}.log"
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
        "annotated/{p}/vcfanno/all.{p}.vep.vcfanno.vcf",
    output:
         db="annotated/{p}/gemini.db",
    log:
        "logs/vcf2db/vcf2db.{p}.log"
    params:
        ped=config["run"]["ped"],
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")


