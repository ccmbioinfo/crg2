from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk FilterMutectCalls "
    "-R {snakemake.input.fasta} -V {snakemake.input.vcf} "
    "--tumor-segmentation {snakemake.input.segments} "
    "--contamination-table {snakemake.input.contamination} "
    "{extra} "
    "-O {snakemake.output.vcf} "
    "{log}"
)