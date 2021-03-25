PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        expand("report/{p}", p=["panel", "panel-flank"]) if config["run"]["panel"] else [],
        "report/coding",
        "report/sv",
        "report/str/{project}_EH_v1.1.xlsx",
        "qc/multiqc/multiqc.html",
        "plots/depths.svg",
        "plots/allele-freqs.svg",
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
