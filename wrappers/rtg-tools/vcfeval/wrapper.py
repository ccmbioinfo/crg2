__author__ = "Aarthi Mohan"
__copyright__ = "Copyright 2021, Aarthi Mohan"
__email__ = "aarthi.mohan@sickkids.ca"
__license__ = "MIT"


import os, pandas as pd

from snakemake.shell import shell

benchmark = pd.read_table(snakemake.input.benchmark, dtype=str).set_index(["pipeline"], drop=False)
pipeline = snakemake.params.get("pipeline", "wes")
baseline = benchmark.loc[pipeline,'baseline_vcf']
eval_bed = benchmark.loc[pipeline,'high_conf_bed']
java_opts = snakemake.params.get("java_opts", "-Xmx10g")
sdf = snakemake.params.sdf

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "export _JAVA_OPTIONS=\"{java_opts}\" &&  "
    "rtg vcfeval "
    "-b {baseline} "
    "-c {snakemake.input.vcf} "
    "-t {sdf} "
    "-e {eval_bed} "
    "-o {snakemake.output} {log} "
)