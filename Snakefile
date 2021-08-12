PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

##### Target rules #####
project = config["run"]["project"]
if config["run"]["pipeline"] == "wes":
    rule all:
        input:
            "report/coding/{}".format(project)
else:
    rule all:
        input:
            expand("report/{p}/{family}", p=["panel", "panel-flank"], family=project) if config["run"]["hpo"] or config["run"]["panel"]  else [],
            expand("report/{p}/{family}", p=["denovo"], family=project) if config["run"]["ped"] else [],
            expand("report/{p}/{family}", p=["panel", "panel-flank", "denovo"], family=project) if (config["run"]["hpo"] or config["run"]["panel"]) and config["run"]["ped"] else [],
            "report/coding/{family}".format(family=project),
            "report/sv",
            "report/str/{family}.EH-v1.1.xlsx".format(family=project),
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

#Common to both wes and wgs
include: "rules/mapping.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/validation.smk"

if config["run"]["pipeline"] == "wes":
    base = "rules/cre/"
else:
    base = "rules/"
    include: base + "sv.smk"
    include: base + "svreport.smk"
    include: base + "str.smk"
    include: base + "validation.smk"

include: base + "calling.smk"
include: base + "filtering.smk"
include: base + "annotation.smk"
include: base + "snvreport.smk"

