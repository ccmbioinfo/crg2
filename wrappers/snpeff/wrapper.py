from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(snpEff {snakemake.params.memory_initial} "
    "{snakemake.params.memory_max} "
    "-i VCF "
    "-o VCF "
    "-dataDir {snakemake.params.data_dir} "
    "{snakemake.params.reference} "
    "{snakemake.input} > {snakemake.output}) {log}"
)
