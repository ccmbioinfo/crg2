__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

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

if not snakemake.threads:
    threads = 8
else:
    threads =  snakemake.threads

shell(
    "bcftools concat {snakemake.params} --threads {threads} {compress_flags} -o {snakemake.output[0]} "
    "{snakemake.input.vcf} {log}"
)
