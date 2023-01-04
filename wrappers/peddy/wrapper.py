from snakemake.shell import shell
import os.path
from os import path


log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

family = snakemake.wildcards.family
dir = os.path.dirname(snakemake.output[0])
prefix = path.join(dir, family)  # qc/peddy/412

shell(
    "peddy {snakemake.params} --prefix {prefix} {snakemake.input.vcf} {snakemake.input.ped} {log} "
)