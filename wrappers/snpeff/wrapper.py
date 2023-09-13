from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# replace symbolic <INS> alleles so that snpeff will annotate against genes
# https://github.com/pcingola/SnpEff/issues/344
replace_ins = "zcat {snakemake.input} | sed 's/<INS>/N/g' | sed 's/<INS:ME:ALU>/N/g' | sed 's/<INS:ME:LINE1>/N/g' | sed 's/<INS:ME:SVA>/N/g' > tmp.vcf; "
gzip_vcf = "bgzip tmp.vcf; "
mv_vcf = "mv tmp.vcf.gz {snakemake.input}; "

shell(
    "(" + replace_ins + gzip_vcf + mv_vcf + "snpEff {snakemake.params.java_opts} "
    "-i VCF "
    "-o VCF "
    "-dataDir {snakemake.params.data_dir} "
    "-s {snakemake.output.report} "
    "{snakemake.params.reference} "
    "{snakemake.input} > {snakemake.output.vcf}) {log}"
)