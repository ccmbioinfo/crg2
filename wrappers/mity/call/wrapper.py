from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
prefix = snakemake.params.prefix
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

mod1 = " module load python/3.9.2_torch_gpu; "
mod2 = " module load freebayes/1.3.1; "
python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity call --prefix {prefix} --out-folder-path {outdir} {snakemake.input.bam} "
try: # hpf
    shell("(" + mod1 + mod2 + python + mity + ") {log}")
except:  # cheo-ri
    shell("(" + python + mity + ") {log}")
