PIPELINE_VERSION="0.9.0"



include: "rules/common.smk"

project=config["run"]["project"]

##### Target rules #####

rule all:
    input:
        "report/{}".format(project)
        
        

##### Modules #####

include: "rules/mapping.smk"
include: "rules/cre/calling.smk"
include: "rules/cre/filtering.smk"
include: "rules/cre/annotation.smk"
include: "rules/cre/snvreport.smk"

# include: "rules/filtering.smk"
# include: "rules/stats.smk"
# include: "rules/qc.smk"
# include: "rules/annotation.smk"
# include: "rules/snvreport.smk"
# include: "rules/sv.smk"
# include: "rules/svreport.smk"
