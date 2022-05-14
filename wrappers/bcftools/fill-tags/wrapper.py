__author__ = "Nour Hanafi"
__copyright__ = "Copyright 2022, Nour Hanafi"
__email__ = "nour.hanafi@sickkids.ca"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


fill_tags_cmd = (
    "bcftools +fill-tags {snakemake.input} " 
    "-o {snakemake.output} -O v "
    "-- -t 'DP2=sum(FORMAT/DP)'; "
)

rename_dp_cmd = (
    "sed -i 's/DP2/DP/' {snakemake.output}"
)

shell("(" + fill_tags_cmd + rename_dp_cmd + ") {log}")
