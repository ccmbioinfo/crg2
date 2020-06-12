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
    threads: 8
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
        known = config["ref"]["known-variants"]
    output:
        bam = protected("recal/{sample}-{unit}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("gatk", "baserecalibrator")

rule remove_decoy:
    input:
        bam = "recal/{sample}-{unit}.bam",
        decoy = config["params"]["remove_decoy"]["bed"],
    output: 
        temp_out = temp("decoy_rm/{sample}-{unit}.no_decoy_reads.temp.bam"),
        out_f = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam"
    log:
        "logs/remove_decoys/{sample}-{unit}.log"
    threads: 8
    resources: mem_mb = 10000
    shell:
        """
        samtools view {input.bam} -b -h -t {threads} -U {output.temp_out} -L {input.decoy}
        samtools view -h -t {threads} {output.temp_out} | grep -v hs37d5 | grep -v NC_007605 | samtools view - -hb > {output.out_f}
        """

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        get_wrapper_path("samtools", "index")
