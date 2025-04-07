from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

python = " export PYTHONPATH={pythonpath}; "

# Temporary fix until mity fixes bug in their code base
# Creating a copy of the input file with a different name to prevent xlsxwriter from throwing 
# "xlsxwriter.exceptions.InvalidWorksheetName: Excel worksheet name '...' must be <= 31 chars." error

orig_input_filename=f"{snakemake.input}"
new_filename=f"mitochondrial_variants/{family}.vcf.gz"
copy_input_file = " cp {orig_input_filename} {new_filename}; "
copy_input_file_tbi = " cp {orig_input_filename}.tbi {new_filename}.tbi; "

shell ("(" + copy_input_file + copy_input_file_tbi + ") {log} ")

mity = " {tool}/mity report --prefix {family} -k --output-dir {outdir} {new_filename};"
#rename_vcf_file= " mv mitochondrial_variants/{family}.mity.normalise.decompose.mity.annotated.vcf mitochondrial_variants/{family}.mity.annotated.vcf ; "

# Remove input file copy
rm_copy= "rm {new_filename}; "
rm_copy_tbi= "rm {new_filename}.tbi;"

shell("(" + python + mity + rm_copy + rm_copy_tbi +") {log}")