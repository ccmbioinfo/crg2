from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
tmpdir = os.path.join(os.getcwd(), snakemake.params.tempdir)


shell(
    "mkdir -p {tmpdir}/kmc"
)

shell(
    "(shovill --outdir {snakemake.params.outdir} --tmpdir {tmpdir} --minlen {snakemake.params.minlength} --R1 {snakemake.input.read1} --R2 {snakemake.input.read2}) {log}"
)