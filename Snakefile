PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"
include: "rules/gather_input.smk"
include: "rules/qc.smk"


sra_run=pd.read_table("accession.tsv", dtype=str)["sra_run"][0]

##### Target rules #####
project = config["run"]["project"]
sample=project

rule all:
    input:
        f"fastq/{sra_run}_1.fastq",
        f"fastq/{sra_run}_2.fastq",
        f"qc/fastqc/{sra_run}",
        f"qc/kraken2/{sra_run}_report.txt",
        f"qc/kraken2/{sra_run}_output.txt",
        f"qc/krona/{sra_run}.taxonomy.krona.html",
        f"qc/bracken/{sra_run}.bracken.txt"