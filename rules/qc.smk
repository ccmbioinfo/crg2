rule fastqc:
    input:
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        outdir="qc/pre-assembly/fastqc/{sra_run}"
    output:
        directory("qc/pre-assembly/fastqc/{sra_run}")
    log:
        "logs/qc/pre-assembly/fastqc/{sra_run}.log"
    wrapper:
        get_wrapper_path("fastqc")


rule kraken2_raw_reads:
    input:
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        kraken2_db = config["tools"]["kraken2_db"],
        qc_type="pre-assembly"
    output:
        report_file="qc/pre-assembly/kraken2_raw_reads/{sra_run}_report-raw_reads.txt",
        output_file="qc/pre-assembly/kraken2_raw_reads/{sra_run}_output-raw_reads.txt"
    log:
        "logs/qc/pre-assembly/kraken2_raw_reads/{sra_run}.log"
    wrapper:
        get_wrapper_path("kraken2")

rule krona_pre_assembly:
    input:
        "qc/pre-assembly/kraken2_raw_reads/{sra_run}_output-raw_reads.txt"
    params:
        qc_type="pre-assembly"
    output:
        "qc/pre-assembly/krona/{sra_run}.taxonomy.krona.pre_assembly.html"
    log:
        "logs/qc/pre-assembly/krona/{sra_run}.log"
    wrapper:
        get_wrapper_path("krona")

rule bracken_pre_assembly:
    input:
        "qc/pre-assembly/kraken2_raw_reads/{sra_run}_report-raw_reads.txt"
    params:
        kraken2_db = config["tools"]["kraken2_db"]
    output:
        "qc/pre-assembly/bracken/{sra_run}.bracken.pre_assembly.txt"
    log:
        "logs/qc/pre-assembly/bracken/{sra_run}.log"
    wrapper:
        get_wrapper_path("bracken")


rule kraken2_post_assembly:
    input:
        assembly_dir="shovill/{sra_run}"
    params:
        kraken2_db = config["tools"]["kraken2_db"],
        qc_type="post-assembly"
    output:
        report_file="qc/post-assembly/kraken2_post_assembly/{sra_run}_report-post_assembly.txt",
        output_file="qc/post-assembly/kraken2_post_assembly/{sra_run}_output-post_assembly.txt"
    log:
        "logs/qc/post-assembly/kraken2_post_assembly/{sra_run}.log"
    wrapper:
        get_wrapper_path("kraken2")

rule krona_post_assembly:
    input:
        "qc/post-assembly/kraken2_post_assembly/{sra_run}_output-post_assembly.txt"
    params:
        qc_type="post-assembly"
    output:
        "qc/post-assembly/krona/{sra_run}.taxonomy.krona.post_assembly.html"
    log:
        "logs/qc/post-assembly/krona/{sra_run}.log"
    wrapper:
        get_wrapper_path("krona")

rule bracken_post_assembly:
    input:
        "qc/post-assembly/kraken2_post_assembly/{sra_run}_report-post_assembly.txt"
    params:
        kraken2_db = config["tools"]["kraken2_db"]
    output:
        "qc/post-assembly/bracken/{sra_run}.bracken.post_assembly.txt"
    log:
        "logs/qc/post-assembly/bracken/{sra_run}.log"
    wrapper:
        get_wrapper_path("bracken")

rule quast:
    input:
        assembly_dir="shovill/{sra_run}",
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        kraken2_db = config["tools"]["kraken2_db"]
    output:
        directory("qc/post-assembly/quast/{sra_run}")
    log:
        "logs/qc/post-assembly/quast/{sra_run}.log"
    wrapper:
        get_wrapper_path("quast")

rule multiqc:
    input:
        [expand(input_file, sra_run=sra_run) for input_file in ["qc/pre-assembly/fastqc/{sra_run}",
                                                "qc/pre-assembly/kraken2_raw_reads/{sra_run}_report-raw_reads.txt",
                                                "qc/pre-assembly/bracken/{sra_run}.bracken.pre_assembly.txt",
                                                "qc/post-assembly/kraken2_post_assembly/{sra_run}_report-post_assembly.txt",
                                                "qc/post-assembly/bracken/{sra_run}.bracken.post_assembly.txt", 
                                                "qc/post-assembly/quast/{sra_run}" 
                                                ]]
    output:
        report = report("qc/multiqc/multiqc.html"),
        stats = "qc/multiqc/multiqc_data/multiqc_general_stats.txt"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")
        
