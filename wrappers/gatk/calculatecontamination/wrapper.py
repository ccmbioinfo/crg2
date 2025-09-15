"""Snakemake wrapper for GATK4 CalculateContamination"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk CalculateContamination  -I {snakemake.input.pileups} "
    "-tumor-segmentation {snakemake.output.segments} "
    "-O {snakemake.output.table} "
    "{log}"
)