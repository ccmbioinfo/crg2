def get_filt_vcf(wildcards):
    if wildcards.p == "gatk_haplotype":
        return "filtered/{family}.gatk_haplotype.vcf.gz"
    else:
        return "filtered/{family}.gatk_mutect2.vcf.gz"


rule vt:
    input: get_filt_vcf
    output:
        "filtered/{p}/{family}.{p}.uniq.normalized.decomposed.vcf"  
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{family}.{p}.log"
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



