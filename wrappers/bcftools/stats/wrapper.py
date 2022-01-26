from snakemake.shell import shell
import os.path
from os import path



extra = snakemake.params.get("extra", "")
vcf = snakemake.input[0]


shell(
    "bcftools stats -s - {vcf} > {snakemake.output[0]}"
)

#rename sample name in output file
dir = os.path.dirname(snakemake.output[0])
base = path.basename(snakemake.output[0]).split('.')[0] 
new_base = base + "_stats.txt"
new_path = path.join(dir, new_base)

shell(
    " sed s/-gatk4//g {snakemake.output[0]} > {new_path} && rm {snakemake.output[0]} "
    "&& mv {new_path} {snakemake.output[0]} ")
