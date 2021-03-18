__author__ = "Aarthi Mohan"
__copyright__ = "Copyright 2021, Aarthi Mohan"
__email__ = "aarthi.mohan@sickkids.ca"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not snakemake.threads:
    threads = 8
else:
    threads =  snakemake.threads

shell(
    "bcftools annotate -a {snakemake.input.annot} -h {snakemake.input.hdr} " 
    "-c CHROM,POS,REF,ALT,INFO/CALLERS,INFO/NUMCALLS "
    "--threads {threads} -o {snakemake.output} -O z "
    "{snakemake.input.vcf} {log}"
)
