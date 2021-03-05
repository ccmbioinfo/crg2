__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from os import path, popen


out_format = "-O v"
outfile = snakemake.output[0]
if outfile.endswith(".vcf.gz"):
    out_format = "-O z"

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.threads:
    threads = snakemake.threads
else:
    threads = 1 

chrom = path.splitext(outfile)[0].split("-")[-1]
chromdir = path.join("samtools", chrom)

region = path.join(snakemake.params.region, chrom + ".txt")    

ref = snakemake.input.ref
bams = snakemake.input.samples
popen("mkdir -p " + chromdir)
command = '"samtools mpileup -t AD -t DP -u -g -f {} {} -r {{}} --BCF --uncompressed | bcftools call -m -v {} > {}/{{}}.vcf "'.format(ref, bams, out_format, chromdir)
#command = '"bcftools mpileup -a AD -a DP -f {} {} -r {{}}   | bcftools call -m -v {} > {}/{{}}.vcf "'.format(ref, bams,out_format, chrom)

shell(
    "cat {region} | parallel -k -j {threads} {command} && "
    "ls {chromdir}/*.vcf > {chrom}.files && "
    "bcftools concat -f {chrom}.files | "
    "bcftools sort > {outfile} && "   
    "rm {chromdir}/*.vcf {chrom}.files && rmdir {chromdir};  {log} "
)