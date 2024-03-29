"""Snakemake wrapper for fastqc."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from os import path
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = path.basename(file_path)

    split_ind = 2 if base.endswith(".gz") else 1
    base = ".".join(base.split(".")[:-split_ind])

    return base


# Run fastqc, since there can be race conditions if multiple jobs
# use the same fastqc dir, we create a temp dir.
with TemporaryDirectory() as tempdir:
    shell(
        "fastqc {snakemake.params} --quiet "
        "--outdir {tempdir} {snakemake.input[0]}"
        " {log}"
    )

    # Move outputs into proper position.
    output_base = basename_without_ext(snakemake.input[0])
    html_path = path.join(tempdir, output_base + "_fastqc.html")
    zip_path = path.join(tempdir, output_base + "_fastqc.zip")
 
   
    if snakemake.output.html != html_path:
        shell("mv {html_path} {snakemake.output.html}")

    if snakemake.output.zip != zip_path:
        shell("mv {zip_path} {snakemake.output.zip}")

    shell(
        " unzip -o {snakemake.output.zip} -d qc/fastqc "
    )

    #rename in qc/fastqc
    unzip_folder = path.join("qc/fastqc", output_base + "_fastqc")
    shell(" mv {unzip_folder} {snakemake.output.unzip_folder}")

    unzip_base = "_".join(output_base.split("_")[:-1])
    fastqc_report_path = path.join(snakemake.output.unzip_folder, "fastqc_data.txt")
    shell("sed -i 's/{output_base}.fastq.gz/{unzip_base}.fastq.gz/g' {fastqc_report_path}")

