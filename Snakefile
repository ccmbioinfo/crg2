PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"

samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)

##### Target rules #####
project = config["run"]["project"]
family=project 

if config["run"]["pipeline"] == "wes":
    rule all:
        input:
            "report/coding/{}".format(project),
            "qc/multiqc/multiqc.html",
            [expand("recal/{family}_{sample}.cram.crai".format(family=project, sample=s)) for s in samples.index],
            [expand("recal/{family}_{sample}.cram.md5".format(family=project, sample=s)) for s in samples.index],
            expand("report/{p}/{family}", p="gatk_somatic", family=project),




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
#  
