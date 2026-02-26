from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "( export PYTHONNOUSERSITE=1 && "
    " defense-finder run -o {snakemake.output} --models-dir {snakemake.params.defensefinder_models} -a {snakemake.input} ) {log}"
)

