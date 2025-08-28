from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(quast.py -o {snakemake.output} -1 {snakemake.input.read1} -2 {snakemake.input.read1} {snakemake.input.assembly_dir}/contigs.fa) {log}"
)