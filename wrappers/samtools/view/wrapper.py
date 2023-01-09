from snakemake.shell import shell


shell(
    "samtools view -T {snakemake.input.ref} -@ {snakemake.threads} -C -o {snakemake.output} {snakemake.input.bam}"
)
