rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}-{unit}.html",
        zip="qc/fastqc/{sample}-{unit}.zip"
    wrapper:
        "file:/hpf/largeprojects/ccm_dccforge/dccdipg/Common/pipelines/crg2/wrappers/fastqc"


rule samtools_stats:
    input:
        "recal/{sample}-{unit}.bam"
    output:
        "qc/samtools-stats/{sample}-{unit}.txt"
    log:
        "logs/samtools-stats/{sample}-{unit}.log"
    wrapper:
        "file:/hpf/largeprojects/ccm_dccforge/dccdipg/Common/pipelines/crg2/wrappers/samtools/stats"


rule multiqc:
    input:
        expand(["qc/samtools-stats/{u.sample}-{u.unit}.txt",
                "qc/fastqc/{u.sample}-{u.unit}.zip",
                "qc/dedup/{u.sample}-{u.unit}.metrics.txt"],
               u=units.itertuples()),
        "snpeff/all.csv"
    output:
        report("qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "logs/multiqc.log"
    wrapper:
        "file:/hpf/largeprojects/ccm_dccforge/dccdipg/Common/pipelines/crg2/wrappers/multiqc"
