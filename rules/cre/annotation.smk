def get_filt_vcf(wildcards):
    if wildcards.p == "coding":
        return "filtered/{family}.vcf.gz"
    elif wildcards.p == "denovo":
        return "filtered/{family}.vcf.gz"
    else:
        return "filtered/{p}/{family}.{p}.vcf.gz".format(p=wildcards.p,family=project)


rule vt:
    input:
        "genotyped/{prefix}.vcf.gz", "genotyped/{prefix}.vcf.gz.tbi"
    output:
        "filtered/{prefix}.uniq.normalized.decomposed.vcf"  
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{prefix}.log"
    wrapper:
        get_wrapper_path("vt")

rule pass:
    input:
       	"{prefix}.{ext}"
    output:
        temp("{prefix}.pass.{ext,(vcf|vcf\.gz)}")
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: 
        filter = "-f PASS"
    wrapper:
        get_wrapper_path("bcftools", "view")

rule vep:
    input:
        "filtered/{p}/{family}.{p}.uniq.normalized.decomposed.pass.vcf",
    output:
        temp("annotated/{p}/vep/{family}.{p}.vep.vcf"),
    log:
        "logs/vep/{family}.vep.{p}.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        ref=config["ref"]["genome"],
    wrapper:
        get_wrapper_path("vep")


rule vcf2db:
    input:
        "annotated/{p}/vep/{family}.{p}.vep.vcf.gz",
    output:
         db="annotated/{p}/{family}-gemini.db",
    log:
        "logs/vcf2db/{family}.vcf2db.{p}.log"
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")


