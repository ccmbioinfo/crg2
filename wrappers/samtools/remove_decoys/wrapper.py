__author__ = "Delvin So"
__copyright__ = "Copyright 2020, Delvin So"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
bam = snakemake.input[0]
decoy = snakemake.input[1]

rm_reads = snakemake.output[0]
out_f = snakemake.output[1]


shell = (
    "samtools view {bam} -b -h -t {snakemake.threads} -o {rm_reads} -U {out_f} -L  {decoy}"
    "samtools view -t {snakemake.threads} {out_f}| grep -v hs37d5 | grep -v NC_007605 | samtools view - -hb > {out_f}"
)