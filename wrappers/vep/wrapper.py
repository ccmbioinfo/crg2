__author__ = "Dennis Kao"
__copyright__ = "Copyright 2020, Dennis Kao"
__email__ = "dennis.kao@sickkids.ca"
__license__ = "BSD"


from snakemake.shell import shell
from os import path
import shutil
import tempfile

shell.executable("bash")

threads = snakemake.threads
threadsprefix = "--fork {}".format(threads) if threads else ""

dir = snakemake.params.get("dir")
dirprefix = "--dir {}".format(dir) if dir else ""

dircache = snakemake.params.get("dir_cache")
dircacheprefix = "--dir_cache {}".format(dircache) if dircache else ""

fasta = snakemake.params.get("ref")
if fasta:
    fastaprefix = "--fasta {}".format(fasta)
else:
    raise Exception("Missing reference fasta file.")

outcalls = snakemake.output[0]
if outcalls.endswith(".vcf.gz"):
    outprefix = "| bcftools view -Oz"
elif outcalls.endswith(".bcf"):
    outprefix = "| bcftools view -Ob"
else:
    outprefix = ""

incalls = snakemake.input[0]
if incalls.endswith(".bcf"):
    incalls = "<(bcftools view {})".format(incalls)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "(vep --vcf -o stdout -i {incalls} {threadsprefix} --species homo_sapiens --no_stats --cache --offline {dirprefix} {dircacheprefix} --symbol --numbers --biotype --total_length "
    "--canonical --gene_phenotype --ccds --uniprot --domains --regulatory --protein --tsl --appris --af --max_af --af_1kg --af_esp --af_gnomad "
    "--pubmed --variant_class --allele_number {fastaprefix} "
    "--plugin SpliceRegion --sift b --polyphen b --hgvs --shift_hgvs 1 --merged  "
    "| sed '/^#/! s/;;/;/g' {outprefix} > {outcalls}) {log}"
)
