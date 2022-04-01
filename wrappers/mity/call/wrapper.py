from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
prefix = snakemake.params.prefix
tool = snakemake.params.tool
print(snakemake.input.bam)

mod1= (" module load python; ")
mod2= (" module load freebayes/1.3.1; ")
mity=(" {tool}/mity call --prefix {prefix} --out-folder-path {outdir} --normalise {snakemake.input.bam} ")
shell ("(" + mod1 + mod2 + mity + ") {log}")