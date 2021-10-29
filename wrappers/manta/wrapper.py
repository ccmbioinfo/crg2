from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(configManta.py "
    "--runDir {snakemake.params.outdir} "
    "--reference {snakemake.input.fasta} "
    "--bam {snakemake.input.bam} "
    "{snakemake.params.include_chroms}; "
    "cd {snakemake.params.outdir}; "
    "./runWorkflow.py "
    "--quiet "
    "-m local "
    "-j {snakemake.threads});"
    "(bcftools filter -i 'FORMAT/SR > 0 & BND_DEPTH >= 10 & MATE_BND_DEPTH >= 10' {snakemake.output.sv} -o {snakemake.output.bnd} -O  z) {log}"
)
