__author__ = "Dennis Kao"
__copyright__ = "Copyright 2020, Dennis Kao"
__email__ = "dennis.kao@sickkids.ca"
__license__ = "BSD"


from snakemake.shell import shell
from os import path, popen
import shutil
import tempfile

shell.executable("bash")

outdb = snakemake.output.db
incalls = snakemake.input[0]

family=snakemake.wildcards.family
ped = snakemake.params.get("ped")
if not ped:
    print('Using a dummy ped for vcf2db.')
    ped = f'{family}.ped'    
    if not path.exists(ped):
       print('Generating a dummy ped for vcf2db.')
       with open(ped, 'w') as f:
          for sample in popen('bcftools query -l {}'.format(incalls)).readlines():
             sample = sample.rstrip() #trailing newline
             f.write('\t'.join(['1', sample, '0', '0', '0', '0']) + '\n')
          
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "(vcf2db.py {incalls} {ped} {outdb}) {log}"
)
