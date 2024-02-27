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

# need sequence dictionary for picard sort
# sometimes Platypus VCFs are not properly sorted , even after vcfstreamsort
seq_dict = snakemake.input.ref.replace(".fa", ".dict")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "platypus callVariants {bams} {params} --nCPU {threads} "
    "--refFile={snakemake.input.ref} {regions} "
    "--output={snakemake.output}  &&  "
    "vcfallelicprimitives -t DECOMPOSED --keep-geno {snakemake.output} "
    "| vcffixup - | vcfstreamsort  > temp && "
    "picard -Xmx2g SortVcf I=temp O=temp.sort.vcf SEQUENCE_DICTIONARY={seq_dict} && "
    "rm temp && mv temp.sort.vcf {snakemake.output} {log} "
)
