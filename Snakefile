PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)

##### Target rules #####
project = config["run"]["project"]
if config["run"]["pipeline"] == "wes":
    rule all:
        input:
            "report/coding/{}".format(project),
            "qc/multiqc/multiqc.html",
            [expand("recal/{family}_{sample}.bam.md5".format(family=config["run"]["project"], sample=s)) for s in samples.index]
elif config["run"]["pipeline"] == "wgs":
    rule all:
        input:
            expand("report/{p}/{family}", p=["panel", "panel-flank"], family=project) if config["run"]["hpo"] or config["run"]["panel"]  else [],
            expand("report/{p}/{family}", p=["denovo"], family=project) if config["run"]["ped"] else [],
            expand("report/{p}/{family}", p=["panel", "panel-flank", "denovo"], family=project) if (config["run"]["hpo"] or config["run"]["panel"]) and config["run"]["ped"] else [],
            "report/coding/{family}".format(family=project),
            "report/sv",
            expand("report/str/{family}.{report_name}.xlsx", family=project, report_name=["EH-v1.1","EHDN"]),
            "qc/multiqc/multiqc.html",
            #"plots/depths.svg",
            #"plots/allele-freqs.svg"
            "programs-{}.txt".format(PIPELINE_VERSION),
            [expand("recal/{family}_{sample}.bam.md5".format(family=config["run"]["project"], sample=s)) for s in samples.index]
elif config["run"]["pipeline"] == "annot":
    rule all:
        input:
            "report/coding/{family}".format(family=project)


localrules: write_version
rule write_version:
    output: "programs-{}.txt".format(PIPELINE_VERSION)
    params: config["tools"]["crg2"]
    shell:
        '''
        sh {params}/get_dependency.sh {params} {output}
        '''

##### Modules #####

if config["run"]["pipeline"] == "wes":
    include: "rules/mapping.smk"
    include: "rules/stats.smk"
    include: "rules/qc.smk"
    base = "rules/cre/"
    include: base + "calling.smk"
    include: base + "filtering.smk"
elif config["run"]["pipeline"] == "wgs":
    include: "rules/mapping.smk"
    include: "rules/stats.smk"
    include: "rules/qc.smk"
    base = "rules/"
    include: base + "calling.smk"
    include: base + "filtering.smk"
    include: base + "sv.smk"
    include: base + "svreport.smk"
    include: base + "str.smk"
elif config["run"]["pipeline"] == "annot":
    base = "rules/"

include: base + "annotation.smk"
include: base + "snvreport.smk"
include: base + "validation.smk"
