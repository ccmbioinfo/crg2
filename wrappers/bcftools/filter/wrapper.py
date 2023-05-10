
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
        filter = soft["freebayes"]["filter"]
        name = soft["freebayes"]["name"]
        command = f"--soft-filter {name} -e '{filter}'  -m + {snakemake.input[0]} "
    elif "samtools" in vcf:
        filter = soft["samtools"]["filter"]
        name = soft["samtools"]["name"]
        command = f"--soft-filter {name} -e '{filter}'  -m + {snakemake.input[0]} "
    elif "platypus" in vcf:
        filter = soft["platypus"]["filter"]
        name = soft["platypus"]["name"]
        command = f"--soft-filter {name} -e '{filter}'  -m + {snakemake.input[0]} "
        command += f" | bcftools filter --soft-filter PASS -i '{filter}' -m + " 
    elif "gatk_haplotype" in vcf:
        snvfilter = soft["gatk"]["snvs"]["filter"]
        snvfilter_name = soft["gatk"]["snvs"]["name"]
        indelfilter = soft["gatk"]["indel"]["filter"]
        indelfilter_name = soft["gatk"]["indel"]["name"]
        command = f"--soft-filter {snvfilter_name} -e '{snvfilter}'  -m + {snakemake.input[0]} "
        command += f" | bcftools filter --soft-filter {indelfilter_name} -e '{indelfilter}'  -m + "
    elif "gatk_somatic" in vcf:
        command = f"{snakemake.input[0]}"
        #shell("mv {snakemake.input[0]} {outfile}")
    else:
        command = ""
else:
    command = ""


shell(
    "bcftools filter {command} {compress_flags} -o {outfile} "

)
