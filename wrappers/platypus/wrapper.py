__author__ = "Aarthi Mohan"
__copyright__ = "Copyright 2020, Aarthi Mohan"
__email__ = "aarthi.mohan@sickkids.ca"
__license__ = "MIT"


import os

from snakemake.shell import shell

regions = snakemake.input.get("regions", "")
if regions:
    regions = "--regions=" + regions

if snakemake.threads:
    threads = snakemake.threads
else:
    threads = 8

bams = snakemake.input.bam
if isinstance(bams, list):
    bams = ",".join(bams)
bams = "--bamFiles={}".format(bams)

params = snakemake.params if snakemake.params else ""

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "platypus callVariants {bams} {params} --nCPU {threads} "
    "--refFile={snakemake.input.ref} {regions} "
    "--output={snakemake.output}  &&  "
    "vcfallelicprimitives -t DECOMPOSED --keep-geno {snakemake.output} "
    "| vcffixup - | vcfstreamsort  > temp && "
    "mv temp {snakemake.output} {log} "
)
