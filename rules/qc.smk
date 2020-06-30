rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}-{unit}.html",
        zip="qc/fastqc/{sample}-{unit}.zip"
    log:
        "logs/fastqc/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("fastqc")


rule samtools_stats:
    input:
        "recal/{sample}-{unit}.bam"
    output:
        "qc/samtools-stats/{sample}-{unit}.txt"
    log:
        "logs/samtools-stats/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("samtools", "stats")


rule fastq_screen:
    input:
        unpack(get_fastq)
    output:
        txt="qc/fastq_screen/{sample}-{unit}_screen.txt",
        png="qc/fastq_screen/{sample}-{unit}_screen.png"
    log:
        "logs/fastq_screen/{sample}-{unit}.log"
    params:
        fastq_screen_config="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/qc/FastQ_Screen_Genomes/fastq_screen.conf", # need to replace bowtie2 bin and path to Arun's directory
        aligner='bowtie2'
    threads: 8
    wrapper:
        get_wrapper_path("fastq_screen")

rule qualimap:
    input: 
        bam = "mapped/{sample}-{unit}.sorted.bam", 
        bai = "mapped/{sample}-{unit}.sorted.bam.bai"
    output:
        "qc/qualimap/{sample}-{unit}/genome_results.txt",
        "qc/qualimap/{sample}-{unit}/qualimapReport.html",
        "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/coverage_histogram.txt"
    params:
        out_dir = "qc/qualimap/{sample}-{unit}",
        nw = config["params"]["qualimap"]["nw"],
        hm = config["params"]["qualimap"]["hm"],
        c = config["params"]["qualimap"]["c"],
        mem_size = config["params"]["qualimap"]["mem"],
        extra = config["params"]["qualimap"]["extra"]
    threads: 8 
    log:
        "logs/qualimap/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("qualimap")


rule verifybamid2:
    input:
        bam = "mapped/{sample}-{unit}.sorted.bam",
        ref = config["ref"]["no_decoy"],
        bai = "mapped/{sample}-{unit}.sorted.bam.bai"
    output:
        sm = "qc/verifybamid2/{sample}-{unit}.selfSM",
    log:
        "logs/verifybamid2/{sample}-{unit}.log"
    params: 
        svd_prefix = config["params"]["verifybamid2"]["svd_prefix"],
        out_dir = "qc/verifybamid2/{sample}-{unit}" ,
        extra = config["params"]["verifybamid2"]["extra"]
    threads: 8
    wrapper:
        get_wrapper_path("verifybamid2")

rule multiqc:
    input:
        expand("qc/samtools-stats/{u.sample}-{u.unit}.txt", u = units.itertuples()),
        expand("qc/fastqc/{u.sample}-{u.unit}.zip", u = units.itertuples()),
        expand("qc/dedup/{u.sample}-{u.unit}.metrics.txt", u = units.itertuples()),
        expand("qc/verifybamid2/{u.sample}-{u.unit}.selfSM", u = units.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/genome_results.txt", u = units.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/coverage_histogram.txt", u = units.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/qualimapReport.html", u = units.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/genome_fraction_coverage.txt", u = units.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt", u = units.itertuples()),
        expand("qc/fastq_screen/{u.sample}-{u.unit}_screen.txt", u = units.itertuples()),
        "snpeff/all.csv"
    params:
    output:
        report("qc/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control") # add dedicated folder
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")


