from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
prefix = os.path.splitext(snakemake.output.json)[0].split(".")

#use module to load EHv5 as it's not yet available via conda
if params.module:
    command = f"module load {params.module} && "
    threads = f" -n {snakemake.get("threads",8)}"
else:
    command = ""
    threads = ""
shell(
    "{command}"
    "sex={snakemake.params.sex} && "
    "( ExpansionHunter"
    " --reference {snakemake.params.ref}"
    " --reads {snakemake.input.bam}"
    " --output-prefix {prefix}"
    " --variant-catalog {snakemake.params.catalog}"
    " --sex $sex 
    " {threads} )"
    " {log}"
)
