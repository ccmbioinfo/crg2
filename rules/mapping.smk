
if config["run"]["input"] == "fastq":
    rule fastq_prep:
        input: 
            units=config["run"]["units"]
        output:
            reads = temp(["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]) 
        log:
            "logs/fastq_prep/{family}_{sample}.log"
        conda:
            "../envs/crg.yaml" 
        script:
            "../scripts/fastq_prep.py"
else:
    rule bamtofastq:
        input:
            bam_file = get_bam
        params:
            outdir = temp("fastq/"),
            sort_check = True
        output:
            fastq1 = temp("fastq/{family}_{sample}_R1.fastq.gz"),
            fastq2 = temp( "fastq/{family}_{sample}_R2.fastq.gz")
        log:
            "logs/bamtofastq/{family}_{sample}.log"
        threads:
            4
        wrapper:
            get_wrapper_path("bedtools", "bamtofastq")



rule map_reads:
    input:
        reads = ["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        temp("mapped/{family}_{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{family}_{sample}.log"
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
        "mapped/{family}_{sample}.sorted.bam"
    output:
        bam = temp("dedup/{family}_{sample}.bam"),
        metrics = "qc/dedup/{family}_{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{family}_{sample}.log"
    params:
        markDuplicates = config["params"]["picard"]["MarkDuplicates"],
        java_opts = config["params"]["picard"]["java_opts"],
    wrapper:
        get_wrapper_path("picard", "markduplicates")


rule recalibrate_base_qualities:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"]
    output:
        bam = protected("recal/{family}_{sample}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/bqsr/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "baserecalibrator")

rule md5:
    input: 
        bam = "recal/{family}_{sample}.bam"
    output:
        md5 = "recal/{family}_{sample}.bam.md5"
    shell:
        """
        md5sum {input.bam} > {output.md5}
        """

rule remove_decoy:
    #redirect first samtools command to a temp output file "decoy.bam" 
    #otherwise it floods the log file with binary stream of decoy reads
    input:
        bam = "recal/{family}_{sample}.bam",
        canon = config["ref"]["canon_bed"],
    output:
        out_f = temp("decoy_rm/{family}_{sample}.no_decoy_reads.bam")
    log:
        "logs/remove_decoys/{family}_{sample}.log"
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

