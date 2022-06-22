__author__ = "Jan Forster"
__copyright__ = "Copyright 2019, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"

from snakemake.shell import shell

## Extract arguments
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if extra == "panel":
    shell(
        " bedtools intersect"
        " -abam {snakemake.input.bams}"
        " -b {snakemake.input.genes}"
        " -ubam > {snakemake.output.bamslice}"
        " && samtools index {snakemake.output.bamslice}"
    )

elif extra == "-header":
    shell(
        "(bedtools intersect"
        " {extra}"
        " -a {snakemake.input.left}"
        " -b {snakemake.input.right}"
        " | bgzip -c > {snakemake.output})"
        " {log}"
    )
