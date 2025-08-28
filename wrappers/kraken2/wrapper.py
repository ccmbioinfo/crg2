from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

print(snakemake.params)

print(snakemake.params.qc_type)

if snakemake.params.qc_type == "post-assembly":
    # Run Kraken2 on post-assembly

    shell(
    "(kraken2 --db {snakemake.params.kraken2_db} --threads 4 --report {snakemake.output.report_file} --output {snakemake.output.output_file} {snakemake.input.assembly_dir}/contigs.fa) {log}"
)

else:
    # Run Kraken2 pre-assembly
    shell(
        "(kraken2 --db {snakemake.params.kraken2_db} --threads 4 --report {snakemake.output.report_file} --output {snakemake.output.output_file} --paired {snakemake.input.read1} {snakemake.input.read2}) {log}"
    )

