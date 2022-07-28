__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell

remove_dups = snakemake.params.markDuplicates
remove_dups_flag = remove_dups.split("=")[1]

if remove_dups_flag == "false":
    # mark duplicates for variant calling
    command = "picard MarkDuplicates {remove_dups} {snakemake.params.java_opts}   INPUT={snakemake.input} OUTPUT={snakemake.output.bam} METRICS_FILE={snakemake.output.metrics} "
else:
    # remove duplicates for the purposes of calculating coverage metrics
    command = "picard MarkDuplicates {remove_dups} {snakemake.params.java_opts} {snakemake.params.validationStringency}  {snakemake.params.assumeSortOrder} INPUT={snakemake.input} OUTPUT={snakemake.output.bam} METRICS_FILE={snakemake.output.metrics} "

print(remove_dups)
shell("(" + command + ") &> {snakemake.log}")
