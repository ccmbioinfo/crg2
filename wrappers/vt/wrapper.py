__author__ = "Dennis Kao"
__copyright__ = "Copyright 2020, Dennis Kao"
__email__ = "dennis.kao@sickkids.ca"
__license__ = "BSD"


from snakemake.shell import shell
from os import path
import shutil
import tempfile

shell.executable("bash")

fasta = snakemake.params.get("ref")
if fasta:
    fastaprefix = "-r {}".format(fasta)
else:
    raise Exception("Missing reference fasta file.")

outcalls = snakemake.output[0]

incalls = snakemake.input[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "(vt decompose -s {incalls} | vt normalize {fastaprefix} -n - | vt uniq - "
    " | vt view -o {outcalls} - ) {log}"
)
