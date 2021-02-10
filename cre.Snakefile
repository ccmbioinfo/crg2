PIPELINE_VERSION="0.9.0"

project=config["run"]["project"]

include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "concat/{}-concat-annot.vcf.gz".format(project)
        

##### Modules #####

include: "rules/mapping.smk"
include: "rules/cre/calling.smk"

# include: "rules/filtering.smk"
# include: "rules/stats.smk"
# include: "rules/qc.smk"
# include: "rules/annotation.smk"
# include: "rules/snvreport.smk"
# include: "rules/sv.smk"
# include: "rules/svreport.smk"
