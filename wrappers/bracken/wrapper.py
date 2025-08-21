from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(bracken -d {snakemake.params.kraken2_db} -i {snakemake.input} -o {snakemake.output} -l S) {log}"
)