from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity call --prefix {family} --output-dir {outdir} -k {snakemake.input.bam} "
shell("(" + python + mity + ") {log}")
