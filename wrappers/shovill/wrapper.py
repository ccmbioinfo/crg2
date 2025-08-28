from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "mkdir -p {snakemake.params.tempdir}"
)

shell(
    "(shovill --outdir {snakemake.params.outdir} --tmpdir {snakemake.params.tempdir} --minlen {snakemake.params.minlength} --R1 {snakemake.input.read1} --R2 {snakemake.input.read2}) {log}"
)