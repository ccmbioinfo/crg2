
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



vcf = snakemake.input[0]
hard = snakemake.params.get("hard", "")
soft = snakemake.params.get("soft", "")
indel = ""

if hard:
    filters = hard
    command = hard + " " + snakemake.input[0]
elif soft:
    if "freebayes" in vcf:
        filter = soft["freebayes"]
        command = "--soft-filter FBQualDepth -e '{}'  -m + {} ".format(filter, snakemake.input[0])
    elif "samtools" in vcf:
        filter = soft["samtools"]
        command = "--soft-filter stQualDepth -e '{}'  -m + {} ".format(filter, snakemake.input[0])
    elif "platypus" in vcf:
        filter = soft["platypus"]
        command = " --soft-filter PlatQualDepth -e '{}'  -m + {} ".format(filter, snakemake.input[0])
    elif "gatk" in vcf:
        snvfilter = soft["gatk"]["snvs"]
        indelfilter = soft["gatk"]["indel"]
        command = "--soft-filter GATKCutoffSNP -e '{}'  -m + {} ".format(snvfilter, snakemake.input[0])
        command += " | bcftools filter --soft-filter GATKCutoffIndel -e '{}'  -m + ".format(indelfilter) 

    else:
        command = ""
else:
    command = ""


shell(
    "bcftools filter {command} {compress_flags} -o {outfile} "

)
