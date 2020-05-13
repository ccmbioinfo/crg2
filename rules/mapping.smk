rule map_reads:
    input:
        reads=get_fastq
    output:
        temp("mapped/{sample}-{unit}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        verbosity=config["params"]["bwa"]["verbosity"],
       	maxMem=config["params"]["bwa"]["maxMem"],
       	markSplitReads=config["params"]["bwa"]["markSplitReads"],
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        get_wrapper_path("bwa", "mem")


rule mark_duplicates:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("dedup/{sample}-{unit}.bam"),
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        get_wrapper_path("picard", "markduplicates")


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"]
    output:
        bam=protected("recal/{sample}-{unit}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
	java_opts=config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("gatk", "baserecalibrator")


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        get_wrapper_path("samtools", "index")
