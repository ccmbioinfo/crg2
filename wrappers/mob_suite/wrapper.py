"""Snakemake wrapper for running mob_suite - a plasmid reconstruction tool."""

__author__ = "Anjali Jain"
__email__ = "anjali.jain@sickkids.ca"


from os import path

from snakemake.shell import shell


assembly = f"{snakemake.input}/contigs.fa"
output_dir = snakemake.output
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(mob_recon --infile {assembly} --outdir {output_dir}) {log}"
)
