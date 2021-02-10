__author__ = "Aarthi Mohan"
__copyright__ = "Copyright 2021, Aarthi Mohan"
__email__ = "aarthi.mohan@sickkids.ca"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not {snakemake.threads}:
    threads =  8
else:
    threads = {snakemake.threads}

shell(
    "bcftools isec -n {snakemake.params.numpass} -O z --threads {threads} -p {snakemake.params.outdir} "
    "{snakemake.input.vcf} {log}"

)
