PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        expand("report/{p}/{family}", p=["panel", "panel-flank"], family=config["run"]["project"]) if config["run"]["panel"] else [],
        "report/coding/{family}".format(family=config["run"]["project"]),
        "report/sv",
        "qc/multiqc/multiqc.html",
        #"plots/depths.svg",
        #"plots/allele-freqs.svg"

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
