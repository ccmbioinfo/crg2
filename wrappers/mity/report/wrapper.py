from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity report --prefix {family} -k --output-dir {outdir} {snakemake.input};"
rename_file= " mv mitochondrial_variants/{family}.mity.normalise.decompose.mity.annotated.vcf mitochondrial_variants/{family}.mity.annotated.vcf"
shell("(" + python + mity + rename_file + ") {log}")