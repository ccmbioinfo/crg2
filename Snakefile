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
            [expand("recal/{family}_{sample}.bam.md5".format(family=config["run"]["project"], sample=s)) for s in samples.index],
            expand("{minio}/{family}", minio=[config["run"]["minio"]], family=project) if config["run"]["minio"] else [],
            expand("report/{p}/{family}", p="gatk_somatic", family=project),
            "report_upload/demultiplexed_reports/{}".format(project) if config["run"]["PT_credentials"] else []
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
            "report/mitochondrial/{family}.mitochondrial.report.csv".format(family=project),
            [expand("recal/{family}_{sample}.bam.md5".format(family=config["run"]["project"], sample=s)) for s in samples.index],
            "report_upload/demultiplexed_reports/{}".format(project) if config["run"]["PT_credentials"] else []

elif config["run"]["pipeline"] == "annot":
    rule all:
        input:
            "report/coding/{family}".format(family=project)
elif config["run"]["pipeline"] == "mity":
    rule all:
        input:
            "report/mitochondrial/{family}.mitochondrial.report.csv".format(family=project)



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
    include: "rules/report_upload.smk"
    base = "rules/cre/"
    include: base + "calling.smk"
    include: base + "filtering.smk"
    include: base + "calling_mosaic.smk"
    include: base + "annotation_mosaic.smk"
    include: base + "snvreport_mosaic.smk"
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
    include: base + "mito_variants.smk"
    include: base + "report_upload.smk"
elif config["run"]["pipeline"] == "annot":
    base = "rules/"
elif config["run"]["pipeline"] == "mity":
    base = "rules/"
    include: base + "mapping.smk"
    include: base + "mito_variants.smk"


include: base + "annotation.smk"
include: base + "snvreport.smk"
include: base + "validation.smk"
