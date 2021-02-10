__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not {snakemake.threads}:
    threads = 8
else:
    threads =  {snakemake.threads}

shell(
    "bcftools concat {snakemake.params} --threads {threads} -o {snakemake.output[0]} "
    "{snakemake.input.vcf} {log}"
)
