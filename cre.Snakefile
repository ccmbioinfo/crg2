PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        expand("called/{project}-{callers}-vt-pass.vcf", project=config["run"]["project"], callers = ["freebayes", "platypus", "samtools_call"]),
        "ensemble/{project}-ensemble.vcf".format(project=config["run"]["project"]),
        #expand("report/{p}", p=["panel", "panel-flank"]) if config["run"]["panel"] else [],
        #"report/coding",
       # "report/sv",
        #"qc/multiqc/multiqc.html",
        #"plots/depths.svg",
        #"plots/allele-freqs.svg"

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
