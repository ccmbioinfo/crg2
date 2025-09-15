"""Snakemake wrapper for GATK4 GetPileupSummaries"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk GetPileupSummaries -I {snakemake.input.cram} "
    "-R {snakemake.input.ref} "
    "-V {snakemake.input.common_variants} "
    "-L {snakemake.input.common_variants} "
    "-O {snakemake.output} "
    "{log}"
)