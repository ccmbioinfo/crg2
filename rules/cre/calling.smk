def get_cre_bams(ext="bam"):
    return ["recal/{}-{}.{}".format(s,units.loc[s].unit[0],ext) if gatk=="gatk" else "recal/gatk3/{}-{}.{}".format(s,units.loc[s].unit[0],ext) for s in samples.index]

rule gatk3:
    input:
        bam=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        known=config["ref"]["known-variants"],
        ref=config["ref"]["genome"],
    output:
        gvcf=protected("called/{project}-gatk3_haplotype.vcf")
    log:
        "logs/gatk/{project}.log"
    params:
        #extra=get_call_variants_params,
        extra=config["params"]["gatk3"]["HaplotypeCaller"],
        annot=config["params"]["gatk3"]["annotation"],
        java_opts=config["params"]["gatk3"]["java_opts"],
    wrapper:
        get_wrapper_path("gatk3", "haplotypecaller")


#duplicating gatk4 rules from crg2/calling.smk for cre file namings 
#sub-workflows does not seem to work smoothly
rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/gatk/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=temp("called/gatk/{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log"
    params:
        extra=get_call_variants_params,
        java_opts=config["params"]["gatk"]["java_opts"],
    group: "gatkcall"
    resources: 
        mem=lambda wildcards, input: len(input.bam) * 15
    wrapper:
        get_wrapper_path("gatk", "haplotypecaller")


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/gatk/{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf=temp("called/gatk/all.{contig}.g.vcf.gz")
    params:
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/combinegvcfs.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "combinegvcfs")


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/gatk/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("called/gatk/genotype/all.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "genotypegvcfs")


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("called/gatk/genotype/all.{contig}.vcf.gz", contig=get_canon_contigs()),
	## use this to remove repetitive contigs for dag generation
	#vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig="GRCh37"), 
    output:
        vcf="called/gatk/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        get_wrapper_path("picard", "mergevcfs")

rule gatk4:
    input: 
        "called/gatk/all.vcf.gz"
    output:
        gvcf=protected("called/{project}-gatk_haplotype.vcf")
    log:
        "logs/gatk/{project}.log"
    shell:
        '''
        gunzip -c -d {input} > {output}
        '''


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
