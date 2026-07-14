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
    resources:
        scratch=20000
    log:
        "logs/shovill/{sra_run}.log"
    wrapper:
        get_wrapper_path("shovill")

rule phager:
    input:
        expand("shovill/{sra_run}/contigs.fa",sra_run=project)
    output:
        directory("phager")
    log:
        expand("logs/phager/{sra_run}.log",sra_run=project)
    wrapper:
        get_wrapper_path("phager")

rule mob_suite:
    input:
        "shovill/{sra_run}/contigs.fa"
    output:
        directory("plasmid_assembly/mob_suite/{sra_run}")
    log:
        "logs/mob_suite/{sra_run}.log"
    wrapper:
        get_wrapper_path("mob_suite")

rule plasmer:
    input:
        "shovill/{sra_run}/contigs.fa"
    output:
        directory("plasmid_assembly/plasmer/{sra_run}")
    params:
        sample_name={sra_run},
        db="/hpf/projects/CTrost/ajain/bacterial_genomics_pipeline/tools/plasmer_db"
    log:
        "logs/plasmer/{sra_run}.log"
    wrapper:
        get_wrapper_path("plasmer")