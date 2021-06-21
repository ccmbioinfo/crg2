PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

##### Target rules #####
import datetime
rule all:
    input:
        expand("report/{p}/{family}", p=["panel", "panel-flank"], family=config["run"]["project"]) if config["run"]["hpo"] or config["run"]["panel"]  else [],
        expand("report/{p}/{family}", p=["denovo"], family=config["run"]["project"]) if config["run"]["ped"] else [],
        expand("report/{p}/{family}", p=["panel", "panel-flank", "denovo"], family=config["run"]["project"]) if (config["run"]["hpo"] or config["run"]["panel"]) and config["run"]["ped"] else [],
        "report/coding/{family}".format(family=config["run"]["project"]),
        "report/sv",
        "report/str/{family}.EH-v1.1.xlsx".format(family=config["run"]["project"]),
        "qc/multiqc/multiqc.html",
        #"plots/depths.svg",
        #"plots/allele-freqs.svg"
        "programs-{}.txt".format(PIPELINE_VERSION)

localrules: write_version
rule write_version:
    output: "programs-{}.txt".format(PIPELINE_VERSION)
    params: config["tools"]["crg2"]
    shell:
        '''
        sh {params}/get_dependency.sh {params} {output}
        '''

##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"
include: "rules/sv.smk"
include: "rules/svreport.smk"
include: "rules/str.smk"
