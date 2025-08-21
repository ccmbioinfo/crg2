from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Take the first two columns of the kraken2 output
shell(
    "cut -f 2,3 {snakemake.input} > qc/krona/{snakemake.wildcards.sra_run}_output.krona"
)

# Run Krona
shell(
    "(ktImportTaxonomy -o {snakemake.output} qc/krona/{snakemake.wildcards.sra_run}_output.krona) {log}"
)