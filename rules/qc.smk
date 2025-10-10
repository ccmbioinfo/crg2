rule fastqc:
    input:
        reads=["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        html="qc/fastqc/{family}_{sample}.html",
        zip="qc/fastqc/{family}_{sample}.zip",
        unzip_folder=directory("qc/fastqc/{family}_{sample}")
    log:
        "logs/fastqc/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("fastqc")


rule samtools_stats:
    input:
        "recal/{family}_{sample}.bam"
    output:
        "qc/samtools-stats/{family}_{sample}.txt"
    log:
        "logs/samtools-stats/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("samtools", "stats")


rule fastq_screen:
    input:
        reads=["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        txt="qc/fastq_screen/{family}_{sample}_screen.txt",
        png="qc/fastq_screen/{family}_{sample}_screen.png"
    log:
        "logs/fastq_screen/{family}_{sample}.log"
    params:
        fastq_screen_config=config["qc"]["fastq_screen"]["conf"], 
        aligner='bowtie2'
    threads: 4
    wrapper:
        get_wrapper_path("fastq_screen")

rule qualimap:
    input: 
        bam = "recal/{family}_{sample}.bam", 
        bai = "recal/{family}_{sample}.bam.bai"
    output:
        "qc/qualimap/{family}_{sample}/genome_results.txt",
        "qc/qualimap/{family}_{sample}/qualimapReport.html",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt"
    params:
        out_dir = "qc/qualimap/{family}_{sample}",
        nw = config["params"]["qualimap"]["nw"],
        hm = config["params"]["qualimap"]["hm"],
        c = config["params"]["qualimap"]["c"],
        mem_size = config["params"]["qualimap"]["mem"],
        pipeline = config["run"]["pipeline"],
        feature = config["params"]["qualimap"]["gtf"],
        extra = config["params"]["qualimap"]["extra"]
    threads: 8
    log:
        "logs/qualimap/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("qualimap")


rule subset:
    input: 
        vcf = get_gatk_vcf
    output:
        sample_vcf = "genotyped/{family}_{sample}-gatk4.vcf.gz"
    params:
        samples = get_sample_order,
        filter = "-f 'PASS,.' "
    log: 
        "logs/bcftools/view/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("bcftools", "view")


rule bcftools_stats:
    input:
        sample_vcf = "genotyped/{family}_{sample}-gatk4.vcf.gz"
    output:
        "qc/bcftools_stats/{family}_{sample}.txt"
    log:
        "logs/bcftools/stats/{family}_{sample}.log"
    params:
    threads: 4
    wrapper:
        get_wrapper_path("bcftools", "stats")


rule multiqc:
    input:
        [expand(input_file, sample=samples.index,family=project) for input_file in ["qc/samtools-stats/{family}_{sample}.txt", 
                                                            "qc/fastqc/{family}_{sample}", 
                                                            "qc/dedup/{family}_{sample}.metrics.txt",  
                                                            "qc/qualimap/{family}_{sample}/genome_results.txt", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt",
                                                            "qc/qualimap/{family}_{sample}/qualimapReport.html", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt", 
                                                            "qc/fastq_screen/{family}_{sample}_screen.txt",
                                                            "qc/bcftools_stats/{family}_{sample}.txt",
                                                                    ]]
    params:
        config["params"]["multiqc"]["config"]
    output:
        report = report("qc/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
        stats = "qc/multiqc/multiqc_data/multiqc_general_stats.txt"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")


rule remove_duplicates:
    # use samtools instead of picard to remove duplicates because picard complains with cram input
    input: 
        cram="recal/{family}_{sample}.cram"
    params:
        ref_cache=config["ref"]["genome"]
    output:
        cram=temp("dups_removed/{family}_{sample}.cram")
    log:
        "logs/remove_duplicates/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("samtools", "markdup")
