from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.params.exclude_chroms:
    excludechroms = "--excludechroms {}".format(snakemake.params.exclude_chroms)
else:
    excludechroms = ""
shell(
    "(export TMPDIR={snakemake.params.outdir}; echo $TMPDIR; "
    "smoove call -x "
    "--name {snakemake.params.name} "
    "--outdir {snakemake.params.outdir} "
    "--fasta {snakemake.input.fasta} "
    "-p {snakemake.threads} "
    "--genotype "
    "{excludechroms} "
    "{snakemake.input.bam}) {log}; "
    "cd {snakemake.params.outdir}; "
    "rm *.bam* *.histo ;"
)
