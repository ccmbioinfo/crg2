rule bamtofastq:
    input:
        bam_file = get_bam
    params:
        outdir = temp("fastq/"),
        sort_check = True
    output:
        fastq1 = temp("fastq/{sample}-{unit}_1.fq"),
        fastq2 = temp("fastq/{sample}-{unit}_2.fq")
    log:
        "logs/bamtofastq/{sample}-{unit}.log"
    threads:
        8
    wrapper:
        get_wrapper_path("bedtools", "bamtofastq")

rule map_reads:
    input:
        reads = get_fastq
    output:
        temp("mapped/{sample}-{unit}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index = config["ref"]["genome"],
        extra = get_read_group,
        verbosity = config["params"]["bwa"]["verbosity"],
        maxMem = config["params"]["bwa"]["maxMem"],
        markSplitReads = config["params"]["bwa"]["markSplitReads"],
        sort = "samtools",
        sort_order = "coordinate"
    threads: 16
    resources: 
        mem=lambda wildcards, threads: threads * 2
    wrapper:
        get_wrapper_path("bwa", "mem")


rule mark_duplicates:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam = temp("dedup/{sample}-{unit}.bam"),
        metrics = "qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        get_wrapper_path("picard", "markduplicates")


rule recalibrate_base_qualities:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"],
        bed = "mapped/{sample}-{unit}-callable.bed" 
    output:
        bam = protected("recal/{sample}-{unit}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("gatk", "baserecalibrator")

rule realignertargetcreator:
    input:
        bam = "mapped/{sample}-{unit}.sorted.bam",
        bai = "mapped/{sample}-{unit}.sorted.bam.bai",
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"]
    output:
        intervals="recal/gatk3/realignertargetcreator/{sample}-{unit}.intervals",
    log:
        "logs/gatk3/realignertargetcreator/{sample}-{unit}.log"
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["RealignerTargetCreator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    threads: 8
    #resources: #cannot access threads here; fix later
     #   mem=lambda wildcards, threads: {threads} * 3
     #using already installed haplotypecaller conda-env, otherwise conda takes forever; 
     #todo: change to a common gatk3.yaml instead of seperate conda for sub-commands
    conda: 
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "realignertargetcreator")

rule indelrealigner:
    input:
        bam = "mapped/{sample}-{unit}.sorted.bam",
        bai = "mapped/{sample}-{unit}.sorted.bam.bai",
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"],
        target_intervals="recal/gatk3/realignertargetcreator/{sample}-{unit}.intervals",
    output:
        bam = protected("recal/gatk3/indelrealigner/{sample}-{unit}-realign.bam"),
    log:
        "logs/gatk3/indelrealigner/{sample}-{unit}.log"
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["IndelRealigner"],
        java_opts = config["params"]["gatk"]["java_opts"],
    conda:
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "indelrealigner")

rule mosdepth:
    input:
        bam = "mapped/{sample}-{unit}.sorted.bam",
        bai = "mapped/{sample}-{unit}.sorted.bam.bai",
    output:
        qbed = "mapped/{sample}-{unit}.quantized.bed.gz",
        bed = "mapped/{sample}-{unit}-callable.bed"
    params:
        prefix = f"mapped/{{sample}}-{{unit}}",
        by = config["ref"]["canon_bed"]
    log:
        "logs/mosdepth/{sample}-{unit}.log"
    shell:
        '''
        export MOSDEPTH_Q0=NO_COVERAGE;
        export MOSDEPTH_Q1=LOW_COVERAGE;
        export MOSDEPTH_Q2=CALLABLE;
        mosdepth -n -F 1804 -Q 1 --quantize 0:1:4: --by {params.by} {params.prefix} {input.bam} 2> {log}
        zgrep "CALLABLE" {output.qbed} > {output.bed}
        '''
       
rule baserecalibrator:
    input:
        bam = "recal/gatk3/indelrealigner/{sample}-{unit}-realign.bam",
        bai = "recal/gatk3/indelrealigner/{sample}-{unit}-realign.bam.bai",
        bed = "mapped/{sample}-{unit}-callable.bed",
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"]
    output:
        "recal/gatk3/baserecalibrator/{sample}-{unit}-recal.grp"
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk3/baserecalibrator/{sample}-{unit}.log"
    threads: 8
    conda:
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "baserecalibrator")

rule printreads:
    input:
        bam = "recal/gatk3/indelrealigner/{sample}-{unit}-realign.bam",
        bai = "recal/gatk3/indelrealigner/{sample}-{unit}-realign.bam.bai",
        ref = config["ref"]["genome"],
        recal_data = "recal/gatk3/baserecalibrator/{sample}-{unit}-recal.grp"
    output:
        protected("recal/gatk3/{sample}-{unit}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["PrintReads"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk3/printreads/{sample}-{unit}.log"
    conda:
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "printreads")

rule remove_decoy:
    #redirect first samtools command to a temp output file "decoy.bam" 
    #otherwise it floods the log file with binary stream of decoy reads
    input:
        bam = "recal/{sample}-{unit}.bam",
        canon = config["ref"]["canon_bed"],
    output:
        out_f = temp("decoy_rm/{sample}-{unit}.no_decoy_reads.bam")
    log:
        "logs/remove_decoys/{sample}-{unit}.log"
    threads: 8
    resources: mem_mb = 10000
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view --threads {threads} -h -L {input.canon} {input.bam} | egrep -v "hs37d5|NC_007605" | samtools view --threads {threads} - -bh > {output.out_f}
        """

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        get_wrapper_path("samtools", "index")
