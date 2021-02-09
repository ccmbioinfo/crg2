def get_cre_bams(ext="bam"):
    return ["recal/{}-{}.{}".format(s,units.loc[s].unit[0],ext) for s in samples.index]

def get_cre_vcfs():
    return ["called/{}-{}-vt.vcf".format(config["run"]["project"],i) for i in ["freebayes", "platypus", "samtools_call"] ]
    

rule freebayes:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        #regions="/path/to/region-file.bed"
    output:
        "called/{project}-freebayes.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/{project}.log"
    params:
        extra=config["params"]["freebayes"],         # optional parameters
        chunksize=100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 8
    wrapper:
        get_wrapper_path("freebayes")

rule platypus:
    input:
	# single or list of bam files
        bam=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        #regions="regions.bed" #remove or empty quotes if not using regions
    output:
	    "called/{project}-platypus.vcf"
    threads: 8
    log:
        "logs/platypus/{project}.log"
    wrapper:
       get_wrapper_path("platypus")

rule samtools_call:
    input:
        samples=get_cre_bams(),
        ref=config["ref"]["genome"],
        bai=get_cre_bams(ext="bam.bai"),
        #regions="regions.bed" #remove or empty quotes if not using regions
    output:
	    "called/{project}-samtools_call.vcf"
    params:
        mpileup = config["params"]["samtools"]["mpileup"],
        call = config["params"]["bcftools"]["call"]
    threads: 8
    log:
        "logs/samtools_call/{project}.log"
    wrapper:
       get_wrapper_path("bcftools", "call")

rule vt:
    input:
        "called/{project}-{caller}.vcf"
    output:
        "called/{project}-{caller}-vt.vcf"  
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{project}-{caller}.log"
    wrapper:
        get_wrapper_path("vt")

rule ensemble:
    input: 
        vcf=get_cre_vcfs(),
        ref=config["ref"]["genome"]
    output:
        "ensemble/{project}-ensemble.vcf"
    threads: 8
    log: "logs/ensemble/{project}.log"
    params: 
        numpass = 2,

    resources:
        mem=lambda wildcards, threads: threads * 2
    wrapper:
        get_wrapper_path("bcbio","variation-recall")
    

rule vcf_isec:
    input:
        vcf =  get_cre_vcfs()
    output:
        outdir = dir("isec")
    params:
        numpass: 1+
    threads: 8
    log: "logs/isec/{project}.log"
    wrapper:
        get_wrapper_path("bcftools","isec")

rule tabix:
    input: "{prefix}.vcf.gz"
    output: "{preifx}.vcf.gz.tbi"
    logs: "logs/tabix/{prefix}.log"
    wrapper:
        get_wrapper_path("tabix")

rule vcf_concat:
    input:
        vcf = [ "isec/000{}.vcf.gz".format(i) for i in range(4) ]
    output:
        "concat/{project}-concat.vcf.gz"
    params: "-d none"
    threads: 8
    wrapper:
        get_wrapper_path("bcftools", "concat")
     
