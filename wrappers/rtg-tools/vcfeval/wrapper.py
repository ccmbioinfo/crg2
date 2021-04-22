__author__ = "Aarthi Mohan"
__copyright__ = "Copyright 2021, Aarthi Mohan"
__email__ = "aarthi.mohan@sickkids.ca"
__license__ = "MIT"


import os

from snakemake.shell import shell

java_opts = snakemake.params.get("java_opts", "-Xmx10g")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "export _JAVA_OPTIONS=\"{java_opts}\" &&  "
    "rtg vcfeval "
    "-b {snakemake.params.baseline} "
    "-c {snakemake.input.vcf} "
    "-t {snakemake.params.sdf} "
    "-e {snakemake.params.eval_bed} "
    "-o {snakemake.output} {log} "
)
