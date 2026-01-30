from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(checkm2 predict -i {snakemake.input} -o {snakemake.output} --database_path {snakemake.params.database_path} --threads 10 ) {log}"
)