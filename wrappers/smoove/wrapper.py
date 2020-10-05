from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "(smoove call -x "
    "--name {snakemake.params.name} "
    "--outdir {snakemake.params.outdir} "
    "--fasta {snakemake.input.fasta} "
    "-p {snakemake.threads} "
    "--genotype "
    "--excludechroms {snakemake.params.exclude_chroms} "
    "{snakemake.input.bam}) {log}; "
    "cd {snakemake.params.outdir}; "
    "rm *.bam* *.histo ;"
)
