include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "filtered/all.vcf.gz",
        "annotated/gemini.db",
        "report",
#        "qc/multiqc.html",
#        "plots/depths.svg",
#        "plots/allele-freqs.svg"

##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"
