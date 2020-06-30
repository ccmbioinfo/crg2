"""Snakemake wrapper for verifybamid2."""

__author__ = "Delvin So"
__copyright__ = "Copyright 2020, Delvin So"
__email__ = "delvin.so@sickkids.ca"
__license__ = "MIT"

import os
from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

out_dir = snakemake.params.get("out_dir")
extra = snakemake.params.get("extra")
svd_prefix = snakemake.params.get("svd_prefix")

bam = snakemake.input[0]
ref = snakemake.input[1]


shell(
    "verifybamid2"
    " {svd_prefix} "
    " --Reference {ref} "
    " --Output {out_dir} "
    " --BamFile {bam} "
    " {extra} "
    "{log}"
)

