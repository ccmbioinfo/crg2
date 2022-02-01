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
        bam = "mapped/{family}_{sample}.sorted.bam", 
        bai = "mapped/{family}_{sample}.sorted.bam.bai"
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


rule verifybamid2:
    input:
        bam = "mapped/{family}_{sample}.sorted.bam",
        ref = config["ref"]["no_decoy"],
        bai = "mapped/{family}_{sample}.sorted.bam.bai"
    output:
        sm = "qc/verifybamid2/{family}_{sample}.selfSM",
    log:
        "logs/verifybamid2/{family}_{sample}.log"
    params: 
        svd_prefix = config["params"]["verifybamid2"]["svd_prefix"],
        out_dir = "qc/verifybamid2/{family}_{sample}" ,
        extra = config["params"]["verifybamid2"]["extra"]
    threads: 4
    wrapper:
        get_wrapper_path("verifybamid2")

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

rule peddy:
    input:
        vcf = get_gatk_vcf,
        ped = dummy_ped
    output:
        "qc/peddy/{family}.html",
        "qc/peddy/{family}.het_check.csv",
        "qc/peddy/{family}.ped_check.csv",
        "qc/peddy/{family}.peddy.ped",
        "qc/peddy/{family}.sex_check.csv",
        "qc/peddy/{family}.background_pca.json"
    params:
        "-p 7"
    log:
        "logs/peddy/{family}.log"
    wrapper:
        get_wrapper_path("peddy")


rule multiqc:
    input:
        [expand(input_file, sample=samples.index,family=project) for input_file in ["qc/samtools-stats/{family}_{sample}.txt", 
                                                            "qc/fastqc/{family}_{sample}", 
                                                            "qc/dedup/{family}_{sample}.metrics.txt", 
                                                            "qc/verifybamid2/{family}_{sample}.selfSM", 
                                                            "qc/qualimap/{family}_{sample}/genome_results.txt", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt",
                                                            "qc/qualimap/{family}_{sample}/qualimapReport.html", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt", 
                                                            "qc/fastq_screen/{family}_{sample}_screen.txt",
                                                            "qc/bcftools_stats/{family}_{sample}.txt",
                                                            "qc/peddy/{family}.html",
                                                            "qc/peddy/{family}.het_check.csv",
                                                            "qc/peddy/{family}.ped_check.csv",
                                                            "qc/peddy/{family}.peddy.ped",
                                                            "qc/peddy/{family}.sex_check.csv",
                                                            "qc/peddy/{family}.background_pca.json"
                                                                    ]]
    params:
    output:
        report("qc/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control") 
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")
