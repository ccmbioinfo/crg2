__author__ = "Dennis Kao"
__copyright__ = "Copyright 2020, Dennis Kao"
__email__ = "dennis.kao@sickkids.ca"
__license__ = "BSD"


from snakemake.shell import shell
from os import path
import shutil
import tempfile

shell.executable("bash")

luascript = snakemake.params.get("lua_script")
if luascript:
    luascriptprefix = "-lua {}".format(luascript)
else:
    luascriptprefix = ""

basepath = snakemake.params.get("base_path")
basepathprefix = "-base-path {}".format(basepath) if basepath else ""

conf = snakemake.params.get("conf")
conf = conf if conf else ""

threads = snakemake.threads
threadsprefix = "-p {}".format(str(threads)) if threads else ""

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
    "(vcfanno {threadsprefix} {luascriptprefix} "
    "{basepathprefix} "
    "{conf} "
    "{incalls} | sed -e 's/Number=A/Number=1/g' {outprefix} > {outcalls}) {log}"
)
