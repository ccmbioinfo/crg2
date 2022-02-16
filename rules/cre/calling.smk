def get_cre_bams(ext="bam"):
    if gatk == "gatk3":
        return expand("recal/gatk3/{family}_{sample}.{ext}", family=project, sample=samples.index, ext=ext)
    return expand("recal/{family}_{sample}.{ext}", family=project, sample=samples.index, ext=ext)


rule gatk3:
    input:
        bam=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        known=config["ref"]["known-variants"],
        ref=config["ref"]["genome"],
        regions="mapped/bed/{family}-sort-callable-{contig}.bed",
        #regions="called/gatk3/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output: "called/gatk3/{family}-{contig}.vcf"
    log:
        "logs/gatk3/{family}-{contig}.log"
    params:
        extra=get_call_variants_params,
        annot=config["params"]["gatk3"]["annotation"],
        java_opts=config["params"]["gatk3"]["java_opts"],
    wrapper:
        get_wrapper_path("gatk3", "haplotypecaller")

rule gathervcf:
    input:
        vcfs = lambda w: expand("called/gatk3/{family}-{contig}.vcf", family=project, contig=get_canon_contigs()),
    output:
        gvcf=protected("genotyped/{family}-gatk3_haplotype.vcf")
    wrapper:
        get_wrapper_path("picard","gathervcfs")


#duplicating gatk4 rules from crg2/calling.smk for cre file namings 
#sub-workflows does not seem to work smoothly
rule call_variants:
    input:
        bam=get_sample_bams,
        #bam=get_cre_bams(),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="mapped/bed/{family}-sort-callable-{contig}.bed",
        #regions="called/gatk/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=temp("called/gatk/{family}_{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{family}_{sample}.{contig}.log"
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
        gvcfs=expand("called/gatk/{{family}}_{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf=temp("called/gatk/{family}.{contig}.g.vcf.gz")
    params:
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/{family}.combinegvcfs.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "combinegvcfs")


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/gatk/{family}.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/gatk/{family}.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
        java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/{family}.genotypegvcfs.{contig}.log"
    group: "gatkcall"
    wrapper:
        get_wrapper_path("gatk", "genotypegvcfs")


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("genotyped/gatk/{{family}}.{contig}.vcf.gz", contig=get_canon_contigs()),
	## use this to remove repetitive contigs for dag generation
	#vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig="GRCh37"), 
    output:
        vcf="genotyped/gatk/{family}.vcf.gz"
    log:
        "logs/picard/{family}-merge-genotyped.log"
    wrapper:
        get_wrapper_path("picard", "mergevcfs")

rule gatk4:
    input: 
        "genotyped/gatk/{family}.vcf.gz".format(family=project)
    output:
        gvcf=protected("genotyped/{family}-gatk_haplotype.vcf")
    log:
        "logs/gatk/{family}.log"
    shell:
        '''
        gunzip -c -d {input} > {output}
        '''


rule freebayes:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        regions="mapped/bed/{family}-sort-callable-{{contig}}.bed".format(family=project)
        #regions=config["ref"]["canon_bed"]
    output:
        temp("called/freebayes/{contig}.vcf")  # either .vcf or .bcf
    log:
        "logs/freebayes/{contig}.log"
    params:
        extra=config["params"]["freebayes"]["call"],         # optional parameters
        chunksize=100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 1
    resources:
        mem=lambda wildcards, threads: threads * 4 if threads > 1 else 20
    wrapper:
        get_wrapper_path("freebayes")

rule merge_freebayes:
    input:
        expand("called/freebayes/{contig}.vcf",contig = get_canon_contigs())
    output:
        protected("genotyped/{family}-freebayes.vcf")
    log:
        "logs/freebayes/{family}-merge.log"
    shell:
        '''
        bcftools concat {input} | bcftools sort > {output}  
        '''
    

rule platypus:
    input:
        bam=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        regions="mapped/{family}-sort-callable.bed"
    output:
	    protected("genotyped/{family}-platypus.vcf")
    params: config["params"]["platypus"]
    threads: 16
    resources:
        mem=lambda wildcards, threads: threads * 4
    log:
        "logs/platypus/{family}.log"
    wrapper:
       get_wrapper_path("platypus")

rule samtools_call:
    input:
        samples=get_cre_bams(),
        ref=config["ref"]["genome"],
        bai=get_cre_bams(ext="bam.bai"),
        region="mapped/bed/{family}-sort-callable-{{contig}}.bed".format(family=project)
    output:
	    temp("called/samtools/{contig}.vcf")
    params:
        mpileup = config["params"]["samtools"]["mpileup"],
        call = config["params"]["bcftools"]["call"],
        contig = lambda w: w.contig
    threads: 8
    resources:
        mem=lambda wildcards, threads: threads * 4
    log:
        "logs/samtools_call/samtools-{contig}.log"
    wrapper:
       get_wrapper_path("bcftools", "call")

rule merge_mpileup:
    input:
        vcfs =lambda w: expand("called/samtools/{contig}.vcf", contig = get_canon_contigs())
        #vcfs =lambda w: expand("called/samtools-{contig}.vcf", contig = "GRCh37d5")
    output:
        protected("genotyped/{family}-samtools.vcf")
    conda:
        "../../envs/common.yaml"
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
    conda:
        "../../envs/common.yaml"

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
    conda:
           "../../envs/common.yaml"
    wrapper:
        get_wrapper_path("tabix")
