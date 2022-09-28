from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
prefix = snakemake.params.prefix
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

mod1 = " module load python; "
mod2 = " module load freebayes/1.3.1; "
python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity call --prefix {prefix} --out-folder-path {outdir} {snakemake.input.bam} "
shell("(" + mod1 + mod2 + python + mity + ") {log}")
