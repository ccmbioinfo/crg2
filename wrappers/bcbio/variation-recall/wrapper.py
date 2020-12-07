__author__ = "Aarthi Mohan"
__copyright__ = "Copyright 2020, Aarthi Mohan"
__email__ = "aarthi.mohan@sickkids.ca"
__license__ = "MIT"


import os

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

callers = ",".join([i.split("-")[1] for i in snakemake.input.vcf])


shell(
    "bcbio-variation-recall ensemble --cores {snakemake.threads} "
    "--numpass {snakemake.params.numpass} "
    "--names {callers} "
    "{snakemake.output} "
    "{snakemake.input.ref} "
    "{snakemake.input.vcf} {log} "
)
