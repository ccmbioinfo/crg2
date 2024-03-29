"""Snakemake wrapper for GATK4 Mutect2"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk Mutect2 "  # Tool and its subprocess
    "--input {snakemake.input.map} "  # Path to input mapping file
    "--output {snakemake.output.vcf} "  # Path to output vcf file
    "--reference {snakemake.input.fasta} "  # Path to reference fasta file
    "--germline-resource {snakemake.input.gnomad_germline} "
    "{log}"  # Logging behaviour
)
