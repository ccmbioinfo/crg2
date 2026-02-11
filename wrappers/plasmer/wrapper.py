"""Snakemake wrapper for running Plasmer - a Plasmid assembly tool."""

__author__ = "Anjali Jain"
__email__ = "anjali.jain@sickkids.ca"


from os import path

from snakemake.shell import shell


assembly = snakemake.input
output_dir = snakemake.output
sample_name = snakemake.params.sample_name
db = snakemake.params.db
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(Plasmer -g {assembly} -p {sample_name} -d {db} -o {output_dir} ) {log}"
)
