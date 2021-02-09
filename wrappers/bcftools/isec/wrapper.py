__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2018, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

if not {snakemake.threads}:
    threads = "--threads " + 8
else:
    threads =  "--threads " + {snakemake.threads}

shell(
    "bcftools isec {snakemake.params.numpass} {threads} -p {snakemake.output.outdir} "
    "{snakemake.input.vcf} "

)
