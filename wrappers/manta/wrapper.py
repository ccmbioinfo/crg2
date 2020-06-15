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
    "-j {snakemake.threads}) {log}"
)
