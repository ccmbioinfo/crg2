__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
params = ""
if {snakemake.params}:
    params = {snakemake.params}


shell("tabix {params} {snakemake.input} {log}")
