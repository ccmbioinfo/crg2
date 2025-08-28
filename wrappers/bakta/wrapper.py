from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(bakta --db {snakemake.params.db} --output {snakemake.params.outdir} --prefix {snakemake.wildcards.sra_run} {snakemake.input}/contigs.fa ) {log}"
)