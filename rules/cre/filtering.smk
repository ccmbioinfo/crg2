def get_cre_vcfs():
    return ["filtered/{}-{}-decomposed-filt.vcf.gz".format(config["run"]["project"],i) for i in [gatk+"_haplotype", "samtools", "freebayes", "platypus"] ]

def get_cre_vcf_tbi():
    return ["filtered/{}-{}-decomposed-filt.vcf.gz.tbi".format(config["run"]["project"],i) for i in [gatk+"_haplotype", "samtools", "freebayes", "platypus"] ]
    
def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="genotyped/all.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{vartype}.log"
    wrapper:
        get_wrapper_path("gatk", "selectvariants")


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="filtered/all.{vartype}.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
        get_wrapper_path("gatk", "variantfiltration")


rule recalibrate_calls:
    input:
        vcf="filtered/all.{vartype}.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.recalibrated.vcf.gz")
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log"
    wrapper:
        get_wrapper_path("gatk", "variantrecalibrator")


rule merge_calls:
    input:
        vcfs=expand("filtered/all.{vartype}.{filtertype}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf="filtered/all.vcf.gz"
    log:
        "logs/picard/merge-filtered.log"
    wrapper:
        get_wrapper_path("picard", "mergevcfs")

rule hard_filter:
    input: 
        "called/{prefix}.vcf.gz", "called/{prefix}.vcf.gz.tbi"
    output: 
        "filtered/{prefix}-hf.vcf.gz"
    params:
        hard = " -i 'ALT=\"<*>\" || QUAL > 5' "
    log:
        "logs/bcftools/hard/{prefix}.log"
    wrapper:
        get_wrapper_path("bcftools", "filter")
        
rule pass:
    input:
        "filtered/{prefix}-hf.vcf.gz", "filtered/{prefix}-hf.vcf.gz.tbi"
    output:
        "filtered/{prefix}-hf-pass.vcf"
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: 
        samples = get_sample_order,
        filter = "-f 'PASS,.' "
    wrapper:
        get_wrapper_path("bcftools", "view")


#rule decompose:
#    input: 
#        "filtered/{prefix}-hf-pass.vcf"
#    output: 
#        "filtered/{prefix}-hf-pass-decomposed.vcf" 
#    log:
#        "logs/vcflib/{prefix}.log"
#    wrapper:
#        get_wrapper_path("vcflib")

rule vt:
    input:
        "filtered/{prefix}-hf-pass.vcf"
    output:
        "filtered/{prefix}-hf-pass-decomposed.vcf"  
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{prefix}.log"
    wrapper:
        get_wrapper_path("vt")

rule soft_filter:
    input: 
        "filtered/{prefix}-hf-pass-decomposed.vcf.gz", "filtered/{prefix}-hf-pass-decomposed.vcf.gz.tbi"
    output:
        "filtered/{prefix}-decomposed-filt.vcf.gz" 
    params: 
        soft = config["filtering"]["soft"],
    log:
        "logs/bcftools/soft/{prefix}.log"
    wrapper:
        get_wrapper_path("bcftools", "filter")


rule vcf_isec:
    input:
        vcf =  get_cre_vcfs(),
        tbi = get_cre_vcf_tbi()
    output:
        vcf = expand("isec/000{index}.vcf.gz", index=range(len(get_cre_vcfs()))),
        sites = "isec/sites.txt"
    params:
        outdir = "isec",
        filters = 'PASS',
        numpass = "+1",
    threads: 8
    log: "logs/isec.log"
    wrapper:
        get_wrapper_path("bcftools","isec")

rule vcf_concat:
    input:
        vcf = expand("isec/000{index}.vcf.gz", index=range(len(get_cre_vcfs())))
    output: 
        "filtered/{project}-decomposed-filt-intersect.vcf.gz"
    params: 
        "-a -d none"
    log: 
        "logs/bcftools/concat/{project}.log"
    threads: 8
    wrapper:
        get_wrapper_path("bcftools", "concat")

rule annot_caller:
    input: "isec/sites.txt"
    output: 
        txt = "isec/sites.caller.txt",
        bz = "isec/sites.caller.txt.gz",
        hdr = "isec/hdr.txt"
    params:
        crg2 = config["tools"]["crg2"]
    
    shell:
        '''
        cat {input} | parallel -k -j 16  {params.crg2}/scripts/annotate-caller.sh {{}} >> {output.txt}
        bgzip -c {output.txt} > {output.bz}
        tabix -s1 -b2 -e2 {output.bz}
        echo -e '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Variant called by">' > {output.hdr}
        '''

rule vcf_annotate:
    input: 
        vcf = "filtered/{project}-decomposed-filt-intersect.vcf.gz",
        annot = "isec/sites.caller.txt.gz",
        hdr = "isec/hdr.txt"
    output: 
        "annotated/{project}-ensemble-decomposed.vcf.gz"
    log: 
        "logs/bcftools/annotate/{project}.log"
    wrapper:
        get_wrapper_path("bcftools", "annotate")



