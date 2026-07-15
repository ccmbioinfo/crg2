from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.params.input_type == "nucleotide":
    input=f"{snakemake.input}/contigs.fa"
else:
    input=f"{snakemake.input}/{snakemake.params.acc}.faa"

shell(
    "( export PYTHONNOUSERSITE=1 && "
    " defense-finder run -o {snakemake.output} --models-dir {snakemake.params.defensefinder_models} -a {input} ) {log}"
)

