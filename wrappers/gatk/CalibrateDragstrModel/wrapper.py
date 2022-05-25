from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk CalibrateDragstrModel "
    "-R {snakemake.input.ref} "
    "-I {snakemake.input.bam} "
    "--str-table-path {snakemake.input.dragstr_table} "
    "-O {snakemake.output[0]} {log}"
)
