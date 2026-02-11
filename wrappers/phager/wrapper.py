"""Snakemake wrapper for running phager - a rapid phage contig predictor using biological feature based machine learning."""

__author__ = "Anjali Jain"
__email__ = "anjali.jain@sickkids.ca"


from snakemake.shell import shell
import sys

assembly = snakemake.input[0]
output_dir = snakemake.output[0]
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
python_path = sys.executable

shell(
    "({python_path} /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/tools/phager-git/phager.py -a {assembly} -o {output_dir} -v) {log}"
)
