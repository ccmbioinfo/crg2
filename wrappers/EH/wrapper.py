from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
prefix = os.path.splitext(snakemake.output.json)[0].split(".")


shell(
    "sex={snakemake.params.sex} && "
    "( ExpansionHunter"
    " --reference {snakemake.params.ref}"
    " --reads {snakemake.input}"
    " --output-prefix {prefix}"
    " --variant-catalog {snakemake.params.catalog}"
    " --sex $sex )"
    " {log}"
)
