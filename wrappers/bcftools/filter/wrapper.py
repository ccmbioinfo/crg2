
from snakemake.shell import shell

outfile = snakemake.output[0]

compress_flags = ""

if outfile.endswith("vcf.gz"):
    compress_flags = "-Oz"
elif outfile.endswith("vcf"):
    compress_flags = "-Ov"
elif outfile.endswith("bcf"):
    compress_flags = "-Ou"
elif outfile.endswith("bcf.gz"):
    compress_flags = "-Ob"


if params.order:
    samples = "-S "
    samples = [samples for sample in popen('bcftools query -l {} | sort'.format(snakemake.input})).readlines()]
    samples = "-S " + ",".join(samples)    
else:
    samples = ""

params = {snakemake.params.soft}
if "freebayes" in {snakemake.input}:
    filter = params["freebayes"]
    filters = "--soft-filter FBQualDepth -e " + filter + " -m + "
elif "samtools" in {snakemake.input}:
    filter = params["samtools"]
    filters = "--soft-filter stQualDepth -e " + filter + " -m + "
elif "platypus" in {snakemake.input}:
    filter = params["platypus"]
    filters = " --soft-filter PlatQualDepth -e " + filter + " -m + "
elif "gatk" in {snakemake.input}:
    snvfilter = params["gatk"]["snv"]
    indelfilter = params["gatk"]["indel"]
    filters = "--soft-filter GATKCutoffSNP -e " + snvfilter + " -m + "
    filters = filters + " --soft-filter GATKCutoffIndel -e" + indelfilter + " -m + "
else:
    filters = ""

shell(
    "bcftools filter {samples} {filters} {snakemake.input[0]} {compress_flags} " "-o {outfile}"
)
