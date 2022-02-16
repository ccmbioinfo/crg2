__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk FilterMutectCalls "
    "-R {snakemake.input.fasta} -V {snakemake.input.vcf} "
    "{extra} "
    "-O {snakemake.output.vcf} "
    "{log}"
)