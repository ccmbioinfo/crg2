

from snakemake.shell import shell




log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")
vcf = snakemake.input[0]


shell(
    "bcftools stats -s - {vcf} > {snakemake.output} {log}"
)