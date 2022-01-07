
from snakemake.shell import shell


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "vcfallelicprimitives  -t DECOMPOSED --keep-geno {snakemake.input[0]} "
    "| vcffixup - | vcfstreamsort "
    " > {snakemake.output[0]}  {log}"
)
