PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"
include: "rules/gather_input.smk"
include: "rules/qc.smk"
include: "rules/assembly.smk"
include: "rules/annotation.smk"


sra_run=pd.read_table("accession.tsv", dtype=str)["sra_run"][0]

##### Target rules #####
project = config["run"]["project"]
sample=project

rule all:
    input:
        f"fastq/{sra_run}_1.fastq",
        f"fastq/{sra_run}_2.fastq",
        f"qc/pre-assembly/fastqc/{sra_run}",
        f"qc/pre-assembly/kraken2_raw_reads/{sra_run}_report-raw_reads.txt",
        f"qc/pre-assembly/kraken2_raw_reads/{sra_run}_output-raw_reads.txt",
        f"qc/pre-assembly/krona/{sra_run}.taxonomy.krona.pre_assembly.html",
        f"qc/pre-assembly/bracken/{sra_run}.bracken.pre_assembly.txt",
        f"shovill/{sra_run}",
        f"qc/post-assembly/kraken2_post_assembly/{sra_run}_report-post_assembly.txt",
        f"qc/post-assembly/kraken2_post_assembly/{sra_run}_output-post_assembly.txt",
        f"qc/post-assembly/quast/{sra_run}",
        f"qc/post-assembly/krona/{sra_run}.taxonomy.krona.post_assembly.html",
        f"qc/post-assembly/bracken/{sra_run}.bracken.post_assembly.txt",
        f"annotated/{sra_run}_bakta",
        f"qc/multiqc/multiqc.html"
