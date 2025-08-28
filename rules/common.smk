
import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version
from datetime import date

min_version("5.7.1")

###### Config file and accession sheets #####
#configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

sra_run=pd.read_table("accession.tsv", dtype=str)["sra_run"][0]
acc_table=pd.read_table("accession.tsv", dtype=str)

validate(acc_table, schema="../schemas/accession.schema.yaml")

project = config["run"]["project"]

##### Wildcard constraints #####
wildcard_constraints:
    sra_run = sra_run,
    sample = str(project)


##### Helper functions #####
    
def get_wrapper_path(*dirs):
    return "file:%s" % os.path.join(workflow.basedir, "wrappers", *dirs)
