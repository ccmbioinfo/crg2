def get_cre_bams(ext="bam"):
    return ["recal/{}-{}.{}".format(s,units.loc[s].unit[0],ext) for s in samples.index]

def get_cre_vcfs():
    return ["called/{}-{}.vcf.gz".format(config["run"]["project"],i) for i in ["freebayes", "platypus", "samtools"] ]

def get_cre_vcf_tbi():
    return ["called/{}-{}.vcf.gz.tbi".format(config["run"]["project"],i) for i in ["freebayes", "platypus", "samtools"] ]
    

rule freebayes:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        #regions="/path/to/region-file.bed"
    output:
        "called/{project}-freebayes.vcf.gz"  # either .vcf or .bcf
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
	    "called/{project}-platypus.vcf.gz"
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
	    "called/{project}-samtools.vcf.gz"
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
        "called/{project}-{caller}.vcf.gz"
    output:
        "called/{project}-{caller}-vt.vcf.gz"  
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
        vcf =  get_cre_vcfs(),
        tbi = get_cre_vcf_tbi()
    output:
        vcf = expand("isec/000{index}.vcf.gz", index=range(len(get_cre_vcfs()))),
        sites = "isec/sites.txt"
    params:
        outdir = "isec",
        numpass = "1+",
    threads: 8
    log: "logs/isec.log"
    wrapper:
        get_wrapper_path("bcftools","isec")

rule tabix:
    input: 
        "{prefix}.vcf.gz"
    output: 
        "{prefix}.vcf.gz.tbi"
    log: 
        "logs/{prefix}.log"
    wrapper:
        get_wrapper_path("tabix")

rule vcf_concat:
    input:
        vcf = expand("isec/000{index}.vcf.gz", index=range(len(get_cre_vcfs())))
    output: 
        "concat/{project}-concat.vcf.gz"
    params: 
        "-a -d none"
    log: 
        "logs/concat.log"
    threads: 8
    wrapper:
        get_wrapper_path("bcftools", "concat")

rule annot_caller:
    input: "isec/sites.txt"
    output: 
        txt = "isec/sites.caller.txt",
        bz = "isec/sites.caller.txt.gz",
        hdr = "isec/hdr.txt"
    shell:
        '''
        cat {input} | parallel -j 16  ./annotate-caller.sh {{}} >> {output.txt}
        bgzip -c {output.txt} > {output.bz}
        tabix -s1 -b2 -e2 {output.bz}
        echo -e '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Variant called by">' > {output.hdr}
        '''

rule vcf_annotate:
    input: 
        vcf = "concat/{project}-concat.vcf.gz",
        annot = "isec/sites.caller.txt.gz",
        hdr = "isec/hdr.txt"
    output: 
        "concat/{project}-concat-annot.vcf.gz"
    log: 
        "logs/annotate.log"
    wrapper:
        get_wrapper_path("bcftools", "annotate")

    


     
