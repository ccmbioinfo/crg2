__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com.com"
__license__ = "MIT"

import os

from snakemake.shell import shell

input_known = snakemake.input.known
extra = snakemake.params.get("extra", "")


bed = snakemake.input.get("bed", None)
if bed is not None:
    bed = "-L " + bed
else:
    bed = ""


input_known_string = ""
if isinstance(input_known, list):
    for known in input_known:
        input_known_string = input_known_string + " -known {}".format(known)
else:
    input_known_string = " -known {}".format(input_known)




output_bai = snakemake.output.get("bai", None)
if output_bai is None:
    extra += " --disable_bam_indexing"


log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "gatk3 {snakemake.params.java_opts} -T IndelRealigner"
    " {extra}"
    " -I {snakemake.input.bam}"
    " -R {snakemake.input.ref}"
    " {input_known_string}"
    " {bed}"
    " --targetIntervals {snakemake.input.target_intervals}"
    " -o {snakemake.output.bam}"
    " {log}"
)
