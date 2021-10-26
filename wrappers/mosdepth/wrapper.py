from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "(export MOSDEPTH_Q0=NO_COVERAGE && "
    " export MOSDEPTH_Q1=LOW_COVERAGE && "
    " export MOSDEPTH_Q2=CALLABLE && "
    " mosdepth -n -F 1804 -Q 1 --quantize {snakemake.params.quantiles} --by {snakemake.params.by} {snakemake.params.prefix} {snakemake.input.bam} && "
    "zgrep 'CALLABLE' {snakemake.output.qbed} > {snakemake.output.bed}) {log}"
)
