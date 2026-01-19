from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


if snakemake.params.user_proteins:
    shell(
        "(/hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/tools/bakta/bin/bakta --db {snakemake.params.db} --output {snakemake.params.outdir} --proteins {snakemake.params.user_proteins} --prefix {snakemake.params.prefix} {snakemake.input}/contigs.fa ) {log}"
    )
# Using the bakta local installation fix to output user proteins into inference.tsv

else:
    shell(
        "(bakta --db {snakemake.params.db} --output {snakemake.params.outdir} --prefix {snakemake.wildcards.sra_run} {snakemake.input}/contigs.fa ) {log}"
    )
