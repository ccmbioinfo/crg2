"""Snakemake wrapper for fastqc."""
__author__ = "Anjali Jain"


from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Run fastqc
shell(
    "mkdir {snakemake.params.outdir}"
)

shell(
        "fastqc --quiet "
        "--outdir {snakemake.params.outdir} {snakemake.input.read1} {snakemake.input.read2}"
        " {log}"
)
