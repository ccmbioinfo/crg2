__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2018, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

input = " ".join(snakemake.input.vcf)

# if one sample, pass this rule, else merge into multisample vcf
if len(snakemake.input.vcf) == 1:
    shell(
    "mv {snakemake.input.vcf} {snakemake.output.vcf_unsort} "
    )

else:
    shell(
    "bcftools merge {extra} -o {snakemake.output.vcf_unsort} "
    "{input} "
    )


shell(
"bcftools sort {snakemake.output.vcf_unsort} "
"-o {snakemake.output.vcf} "
)
