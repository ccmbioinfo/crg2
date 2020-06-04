from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "smoove call -x "
    "--name {snakemake.params.name} "
    #" --exclude $bed "
    "--fasta {snakemake.input.fasta} "
    "-p {snakemake.threads} "
    "--genotype {snakemake.input.bam}"
)
