from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(whamg "
    "-x {snakemake.threads} "
    "-a {snakemake.input.fasta} "
    "-f {snakemake.input.bam} "
    "-c {snakemake.params.include_chroms} > {snakemake.output}) {log} "
)
