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


chrom = snakemake.params.contig

#use callable bed file instead of whole genome
region = snakemake.input.region


ref = snakemake.input.ref
bams = snakemake.input.samples


##command string for using -r
#command = '"samtools mpileup -t AD -t DP -u -g -f {} {} -r {{}} --BCF --uncompressed | bcftools call -m -v {} > {}/{{}}.vcf "'.format(ref, bams, out_format, chromdir)
##command string for bcftools mpileup with -r
#command = '"bcftools mpileup -a AD -a DP -f {} {} -r {{}}   | bcftools call -m -v {} > {}/{{}}.vcf "'.format(ref, bams,out_format, chrom)
##command string for -l bed option (tried parallelizing over each interval in bed, which creates ~100K files causing process to quit!)
command = "samtools mpileup -t AD -t DP -u -g -f {} {} -l {} --BCF --uncompressed | bcftools call -m -v {} | bcftools sort {} -o {}".format(ref, bams,region, out_format, out_format, outfile)
    
shell ("{command}; {log}")