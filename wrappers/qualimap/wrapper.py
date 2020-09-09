"""Snakemake wrapper for qualimap."""

__author__ = "Delvin So"
__copyright__ = "Copyright 2020, Delvin So"
__email__ = "delvin.so@sickkids.ca"
__license__ = "MIT"

import os
from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

num_windows = snakemake.params.get("nw", "")
min_homopolymer = snakemake.params.get("hm", "")
mem_size = snakemake.params.get("mem_size", "")
c = snakemake.params.get("c", "")
extra = snakemake.params.get("extra")

out_dir = snakemake.params.get("out_dir")


shell(
    "qualimap bamqc"
    " -bam {snakemake.input}"
    "  {c}"
    " -nw {num_windows}"
    " -hm {min_homopolymer}"
    " -nt {snakemake.threads}"
    " --java-mem-size={mem_size}"
    " -outdir {out_dir}"
    " {extra}"
    "{log}"
)



