__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"


import os

from snakemake.shell import shell

known = snakemake.input.get("known", "")
if known:
    known = "--dbsnp " + known

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")
bams = snakemake.input.bam
if isinstance(bams, str):
    bams = [bams]
bams = list(map("-I {}".format, bams))

annot = snakemake.params.get("annot","")
if annot:
    extra += " ".join([ " --annotation " + i for i in annot.split(" ") ])

#print(extra)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk3 {java_opts} -T HaplotypeCaller {extra} "
    "-R {snakemake.input.ref} "
    "{bams} "
    "{known} "
    "-o {snakemake.output} {log}"
)
