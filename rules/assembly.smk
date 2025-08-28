rule shovill:
    input:
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        minlength="1000",
        #kmers="2",
        tempdir="temp",
        outdir="shovill/{sra_run}"
    output:
        directory("shovill/{sra_run}")
    log:
        "logs/shovill/{sra_run}.log"
    wrapper:
        get_wrapper_path("shovill")