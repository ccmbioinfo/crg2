def get_cre_bams(ext="bam"):
    return ["recal/{}-{}.{}".format(s,units.loc[s].unit[0],ext) for s in samples.index]

def get_cre_vcfs():
    return ["filtered/{}-{}-vt.vcf.gz".format(config["run"]["project"],i) for i in ["gatk_haplotype", "samtools", "freebayes", "platypus"] ]

def get_cre_vcf_tbi():
    return ["filtered/{}-{}-vt.vcf.gz.tbi".format(config["run"]["project"],i) for i in ["gatk_haplotype", "samtools", "freebayes", "platypus"] ]
    
rule gatk3:
    input:
        bam=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        known=config["ref"]["known-variants"],
        ref=config["ref"]["genome"],
    output:
        gvcf=protected("called/{project}-gatk_haplotype.vcf")
    log:
        "logs/gatk/{project}.log"
    params:
        #extra=get_call_variants_params,
        extra=config["params"]["gatk3"]["HaplotypeCaller"],
        annot=config["params"]["gatk3"]["annotation"],
        java_opts=config["params"]["gatk3"]["java_opts"],
    wrapper:
        get_wrapper_path("gatk3", "haplotypecaller")

        

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
    threads: 16
    resources:
        mem=lambda wildcards, threads: threads * 4
    wrapper:
        get_wrapper_path("freebayes")

rule platypus:
    input:
        bam=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        #regions="regions.bed" #remove or empty quotes if not using regions
    output:
	    "called/{project}-platypus.vcf"
    params: config["params"]["platypus"]
    threads: 16
    resources:
        mem=lambda wildcards, threads: threads * 4
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
	    "called/samtools-{contig}.vcf"
    params:
        mpileup = config["params"]["samtools"]["mpileup"],
        call = config["params"]["bcftools"]["call"],
        region = config["ref"]["split_genome"]
    threads: 8
    resources:
        mem=lambda wildcards, threads: threads * 4
    log:
        "logs/samtools_call/samtools-{contig}.log"
    wrapper:
       get_wrapper_path("bcftools", "call")

rule merge_mpileup:
    input:
        vcfs =lambda w: expand("called/samtools-{contig}.vcf", contig = get_canon_contigs())
        #vcfs =lambda w: expand("called/samtools-{contig}.vcf", contig = "GRCh37d5")
    output:
        "called/{project}-samtools.vcf"
    shell:
        '''
        if [ -f files ]; then rm files; fi;
        for i in {input}; do echo $i >> files; done;
        bcftools concat -f files | bcftools sort > {output}
        rm {input}
        '''
 
rule bgzip:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    shell:
        '''
        bgzip -c {input} > {output}
        '''

rule tabix:
    input: 
        "{prefix}.vcf.gz"
    output: 
        "{prefix}.vcf.gz.tbi"
    log: 
        "logs/{prefix}.log"
    wrapper:
        get_wrapper_path("tabix")
