__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from os import popen

outfile = snakemake.output[0]

compress_flags = ""

if outfile.endswith("vcf.gz"):
    compress_flags = "-Oz"
elif outfile.endswith("vcf"):
    compress_flags = "-Ov"
elif outfile.endswith("bcf"):
    compress_flags = "-Ou"
elif outfile.endswith("bcf.gz"):
    compress_flags = "-Ob"

family = snakemake.wildcards.family
sample = snakemake.wildcards.sample

if snakemake.params.get("samples"): 
    samples = f" -s {family}_{sample} ".format(family=family, sample = sample) 
else: 
    samples = "" 

params = snakemake.params.filter if snakemake.params.get("filter") else ""


shell(
    "bcftools view {samples} {params} {snakemake.input[0]} {compress_flags} -o {outfile} "
)
