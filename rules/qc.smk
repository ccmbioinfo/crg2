rule fastqc:
    input:
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        outdir="qc/fastqc/{sra_run}"
    output:
        directory("qc/fastqc/{sra_run}")
    log:
        "logs/qc/fastqc/{sra_run}.log"
    wrapper:
        get_wrapper_path("fastqc")


rule kraken2:
    input:
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        kraken2_db = config["tools"]["kraken2_db"]
    output:
        report_file="qc/kraken2/{sra_run}_report.txt",
        output_file="qc/kraken2/{sra_run}_output.txt"
    log:
        "logs/qc/kraken2/{sra_run}.log"
    wrapper:
        get_wrapper_path("kraken2")

rule krona:
    input:
        "qc/kraken2/{sra_run}_output.txt"
    output:
        "qc/krona/{sra_run}.taxonomy.krona.html"
    log:
        "logs/qc/krona/{sra_run}.log"
    wrapper:
        get_wrapper_path("krona")

rule bracken:
    input:
        "qc/kraken2/{sra_run}_report.txt"
    params:
        kraken2_db = config["tools"]["kraken2_db"]
    output:
        "qc/bracken/{sra_run}.bracken.txt"
    log:
        "logs/qc/bracken/{sra_run}.log"
    wrapper:
        get_wrapper_path("bracken")


    