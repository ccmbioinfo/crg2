from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

print(snakemake.input.read1)
print(snakemake.input.read2)
print(snakemake.params.kraken2_db)

shell(
    "(kraken2 --db {snakemake.params.kraken2_db} --threads 4 --report {snakemake.output.report_file} --output {snakemake.output.output_file} --paired {snakemake.input.read1} {snakemake.input.read2}) {log}"
)