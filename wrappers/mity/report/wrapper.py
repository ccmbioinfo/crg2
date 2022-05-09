from snakemake.shell import shell
import os
import process_report

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
prefix = snakemake.params.prefix
tool = snakemake.params.tool

mod1 = " module load python; "
mod2 = " module load freebayes/1.3.1; "
bgzip = " bgzip {snakemake.input}; "
mity = " {tool}/mity report --prefix {prefix} --out-folder-path {outdir} {snakemake.input}.gz; "
remove_excel = " rm report/mitochondrial/{family}_mito.annotated_variants.xlsx; "
shell("(" + mod1 + mod2 + bgzip + mity + remove_excel + ") {log}")

report = f"report/mitochondrial/{family}_mito.annotated_variants.csv"
process_report.main(report, family)

remove_csv = " rm report/mitochondrial/{family}_mito.annotated_variants.csv; "
shell("(" + remove_csv + ") {log}")
