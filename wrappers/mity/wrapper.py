from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
prefix = snakemake.params.prefix
tool = snakemake.params.tool
vcfanno_config= snakemake.params.vcfanno_config
pythonpath = tool.replace("bin", "")

python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity runall --prefix {prefix} --output-dir {outdir} -k {snakemake.input.bam}; "
shell("(" + python + mity + ") {log}")

bgzip = " bgzip mitochondrial_variants/{family}.normalise.mity.annotated.vcf;"
tbi = "tabix mitochondrial_variants/{family}.normalise.mity.annotated.vcf.gz "
shell("(" + bgzip + tbi + ") {log}")