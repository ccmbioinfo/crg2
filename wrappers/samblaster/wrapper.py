from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

splitters = snakemake.output.splitters.replace('bam', 'sam')
discordants = snakemake.output.discordants.replace('bam', 'sam')
shell(
    "samtools view -h {snakemake.input} "
    "| samblaster "
    "--acceptDupMarks "
    "--excludeDups "
    "--addMateTags -M "
    "--splitterFile splitters "
    "--discordantFile discordants "
    "--ignoreUnmated "
    "-o /dev/null; "
    "samtools view -Sb splitters > {snakemake.output.splitters}; "
    "samtools view -Sb discordants > {snakemake.output.discordants}"
)
