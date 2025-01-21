from snakemake.shell import shell
import os
import process_report

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
prefix = snakemake.params.prefix
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")

mod1 = " module load python/3.9.2_torch_gpu; "
mod2 = " module load freebayes/1.3.1; "
python = " export PYTHONPATH={pythonpath}; "
bgzip = " bgzip {snakemake.input}; "
mity = " {tool}/mity report --prefix {prefix} --out-folder-path {outdir} {snakemake.input}.gz; "
remove_excel = " rm report/mitochondrial/{family}_mito.annotated_variants.xlsx; "

try: # hpf
    shell("(" + mod1 + mod2 + python + mity + remove_excel + ") {log}")
except:
    shell("(" + python + bgzip + mity + remove_excel + ") {log}")

report = f"report/mitochondrial/{family}_mito.annotated_variants.csv"
process_report.main(f"{snakemake.input}.gz", report, family)

remove_csv = " rm report/mitochondrial/{family}_mito.annotated_variants.csv; "
shell("(" + remove_csv + ") {log}")
