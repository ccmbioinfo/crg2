rule fastqc:
    input:
        read1="fastq/{sra_run}_1.fastq",
        read2="fastq/{sra_run}_2.fastq"
    params:
        outdir="qc/fastqc/{sra_run}"
    output:
        directory("qc/fastqc/{sra_run}")
    log:
        "logs/qc/fastqc/{sra_run}.log"
    wrapper:
        get_wrapper_path("fastqc")
''''
rule multiqc:
    input:
        [expand(input_file, sample=samples.index,family=project) for input_file in ["qc/samtools-stats/{family}_{sample}.txt", 
                                                            "qc/fastqc/{family}_{sample}", 
                                                            "qc/dedup/{family}_{sample}.metrics.txt", 
                                                            "qc/verifybamid2/{family}_{sample}.selfSM", 
                                                            "qc/qualimap/{family}_{sample}/genome_results.txt", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt",
                                                            "qc/qualimap/{family}_{sample}/qualimapReport.html", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt", 
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
                                                            "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt", 
                                                            "qc/fastq_screen/{family}_{sample}_screen.txt",
                                                            "qc/bcftools_stats/{family}_{sample}.txt",
                                                            "qc/peddy/{family}.html",
                                                            "qc/peddy/{family}.vs.html",
                                                            "qc/peddy/{family}.het_check.csv",
                                                            "qc/peddy/{family}.ped_check.csv",
                                                            "qc/peddy/{family}.peddy.ped",
                                                            "qc/peddy/{family}.sex_check.csv",
                                                            "qc/peddy/{family}.background_pca.json"
                                                                    ]]
    params:
        config["params"]["multiqc"]["config"]
    output:
        report = report("qc/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
        stats = "qc/multiqc/multiqc_data/multiqc_general_stats.txt"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")

'''
