from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

family = snakemake.wildcards.family
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

mod1 = " module load python; "
mod2 = " module load freebayes/1.3.1; "
python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity normalise --outfile {snakemake.output} {snakemake.input}"
shell("(" + mod1 + mod2 + python + mity + ") {log}")
