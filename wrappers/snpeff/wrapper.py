from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(snpEff {snakemake.params.java_opts} "
    "-i VCF "
    "-o VCF "
    "-dataDir {snakemake.params.data_dir} "
    "-s {snakemake.output.report} "
    "{snakemake.params.reference} "
    "{snakemake.input} > {snakemake.output.vcf}) {log}"
)
