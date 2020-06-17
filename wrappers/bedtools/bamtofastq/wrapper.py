from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

outdir = snakemake.output.fastq1.split('/')[0]

shell(

    "(samtools sort -n  "
    "-@ {snakemake.threads} "
    "-T {outdir}/{snakemake.wildcards.sample} "
    "{snakemake.input} "
    "| bedtools bamtofastq -i /dev/stdin "
    "-fq {snakemake.output.fastq1} "
    "-fq2 {snakemake.output.fastq2}) {log}"
)
