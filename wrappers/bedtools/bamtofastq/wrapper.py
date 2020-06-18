from snakemake.shell import shell
import pysam

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

sort_cmd = ""

bamtofastq_cmd = (
    " bedtools bamtofastq -i {snakemake.input.bam_file} "
    "-fq {snakemake.output.fastq1} "
    "-fq2 {snakemake.output.fastq2}"
)

# if bam is not sorted by query name, it must be re-sorted
if snakemake.params.sort_check:
    bam = pysam.AlignmentFile(snakemake.input.bam_file, "rb")
    header = bam.header
    sort_order = header["HD"]["SO"]
    bam.close()

    if sort_order != "queryname":
        sort_cmd = (
            "samtools sort -n "
            "-@ {snakemake.threads} "
            "-T {snakemake.params.outdir}/{snakemake.wildcards.sample} "
            "{snakemake.input} |"
        )
        bamtofastq_cmd = (
            " bedtools bamtofastq -i /dev/stdin "
            "-fq {snakemake.output.fastq1} "
            "-fq2 {snakemake.output.fastq2}"
        )

shell(
    "(" + sort_cmd + bamtofastq_cmd + ") {log}"
)
