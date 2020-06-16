from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "(run_metasv.py "
    "--sample {snakemake.params.sample} "
    "--reference {snakemake.input.fasta} "
    "--bam {snakemake.input.bam} "
    "--outdir {snakemake.params.outdir} "
    "--lumpy_vcf {snakemake.input.lumpy} "
    "--manta_vcf {snakemake.input.manta} "
    "--wham_vcf {snakemake.input.wham} "
    "--num_threads {snakemake.threads} "
    "--disable_assembly "
    "--filter_gaps) {log}"
)
