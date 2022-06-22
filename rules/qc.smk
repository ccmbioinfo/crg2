rule fastqc:
    input:
        reads=["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        html="qc/fastqc/{family}_{sample}.html",
        zip="qc/fastqc/{family}_{sample}.zip",
        unzip_folder=directory("qc/fastqc/{family}_{sample}")
    log:
        "logs/fastqc/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("fastqc")


rule samtools_stats:
    input:
        "recal/{family}_{sample}.bam"
    output:
        "qc/samtools-stats/{family}_{sample}.txt"
    log:
        "logs/samtools-stats/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("samtools", "stats")


rule fastq_screen:
    input:
        reads=["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        txt="qc/fastq_screen/{family}_{sample}_screen.txt",
        png="qc/fastq_screen/{family}_{sample}_screen.png"
    log:
        "logs/fastq_screen/{family}_{sample}.log"
    params:
        fastq_screen_config=config["qc"]["fastq_screen"]["conf"], 
        aligner='bowtie2'
    threads: 4
    wrapper:
        get_wrapper_path("fastq_screen")

rule qualimap:
    input: 
        bam = "mapped/{family}_{sample}.sorted.bam", 
        bai = "mapped/{family}_{sample}.sorted.bam.bai"
    output:
        "qc/qualimap/{family}_{sample}/genome_results.txt",
        "qc/qualimap/{family}_{sample}/qualimapReport.html",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt"
    params:
        out_dir = "qc/qualimap/{family}_{sample}",
        nw = config["params"]["qualimap"]["nw"],
        hm = config["params"]["qualimap"]["hm"],
        c = config["params"]["qualimap"]["c"],
        mem_size = config["params"]["qualimap"]["mem"],
        pipeline = config["run"]["pipeline"],
        feature = config["params"]["qualimap"]["gtf"],
        extra = config["params"]["qualimap"]["extra"]
    threads: 8
    log:
        "logs/qualimap/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("qualimap")


rule verifybamid2:
    input:
        bam = "mapped/{family}_{sample}.sorted.bam",
        ref = config["ref"]["no_decoy"],
        bai = "mapped/{family}_{sample}.sorted.bam.bai"
    output:
        sm = "qc/verifybamid2/{family}_{sample}.selfSM",
    log:
        "logs/verifybamid2/{family}_{sample}.log"
    params: 
        svd_prefix = config["params"]["verifybamid2"]["svd_prefix"],
        out_dir = "qc/verifybamid2/{family}_{sample}" ,
        extra = config["params"]["verifybamid2"]["extra"]
    threads: 4
    wrapper:
        get_wrapper_path("verifybamid2")

rule subset:
    input: 
        vcf = get_gatk_vcf
    output:
        sample_vcf = "genotyped/{family}_{sample}-gatk4.vcf.gz"
    params:
        samples = get_sample_order,
        filter = "-f 'PASS,.' "
    log: 
        "logs/bcftools/view/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("bcftools", "view")


rule bcftools_stats:
    input:
        sample_vcf = "genotyped/{family}_{sample}-gatk4.vcf.gz"
    output:
        "qc/bcftools_stats/{family}_{sample}.txt"
    log:
        "logs/bcftools/stats/{family}_{sample}.log"
    params:
    threads: 4
    wrapper:
        get_wrapper_path("bcftools", "stats")

rule peddy:
    input:
        vcf = get_gatk_vcf,
        ped = peddy_ped
    output:
        "qc/peddy/{family}.html",
        "qc/peddy/{family}.vs.html",
        "qc/peddy/{family}.het_check.csv",
        "qc/peddy/{family}.ped_check.csv",
        "qc/peddy/{family}.peddy.ped",
        "qc/peddy/{family}.sex_check.csv",
        "qc/peddy/{family}.background_pca.json"
    params:
        "-p 7"
    log:
        "logs/peddy/{family}.log"
    wrapper:
        get_wrapper_path("peddy")


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


rule remove_duplicates:
    input: 
        "recal/{family}_{sample}.bam"
    params:
        markDuplicates="REMOVE_DUPLICATES=true",
        java_opts = config["params"]["picard"]["java_opts"],
        validationStringency=config["params"]["picard"]["ValidationStringency"],
        assumeSortOrder=config["params"]["picard"]["AssumeSortOrder"]
    output:
        bam="dups_removed/{family}_{sample}.bam",
        metrics="dups_removed/{family}_{sample}.dup.metrics.txt"
    log:
        "logs/remove_duplicates/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("picard", "markduplicates")


rule calculate_coverage:
    input:
        bam="dups_removed/{family}_{sample}.bam",
        bed=config["annotation"]["svreport"]["exon_bed"],
        bed_index=config["ref"]["bed_index"]
    params:
        crg2_dir=config['tools']['crg2']
    output:
        coverage_dir=directory("coverage/{family}_{sample}/"),
    log:
        "logs/coverage/{family}_{sample}.log"
    conda:
        "../envs/coverage.yaml"
    shell:
        """
        mkdir -p coverage/{wildcards.family}_{wildcards.sample}/
        cd coverage/{wildcards.family}_{wildcards.sample}/

        bedtools coverage -d -sorted -a {input.bed} -b ../../{input.bam} -g {input.bed_index} > {wildcards.family}_{wildcards.sample}.dedup.cov

        python {params.crg2_dir}/utils/bam.coverage.base_wise.stat.py {wildcards.family}_{wildcards.sample}.dedup.cov 5 $'\\t' > {wildcards.family}_{wildcards.sample}.coverage_stats.csv

        echo {wildcards.family}_{wildcards.sample}.dedup.cov","`cat {wildcards.family}_{wildcards.sample}.dedup.cov | awk -F '\\t' '{{if ($6<20) {{print $4","$1":"$2+$5","$6}}}}' | wc -l` > {wildcards.family}_{wildcards.sample}.less20x.stats.csv
        {params.crg2_dir}/utils/bam.coverage.less20.sh {wildcards.family}_{wildcards.sample}.dedup.cov > {wildcards.family}_{wildcards.sample}.less20x_coverage.csv

        {params.crg2_dir}/utils/bam.coverage.per_exon.sh {wildcards.family}_{wildcards.sample}.dedup.cov > {wildcards.family}_{wildcards.sample}.per_exon_coverage.csv
        cat {wildcards.family}_{wildcards.sample}.per_exon_coverage.csv | sed 1d  > {wildcards.family}_{wildcards.sample}.per_exon_coverage.csv.tmp
        python {params.crg2_dir}/utils/bam.coverage.base_wise.stat.py {wildcards.family}_{wildcards.sample}.per_exon_coverage.csv.tmp 2 ',' > {wildcards.family}_{wildcards.sample}.per_exon.distribution.csv
        rm {wildcards.family}_{wildcards.sample}.per_exon_coverage.csv.tmp

        median_line=`cat {wildcards.family}_{wildcards.sample}.dedup.cov | wc -l`
        median_line=$(($median_line/2))
        cat {wildcards.family}_{wildcards.sample}.dedup.cov | awk '{{print $6}}' | sort -n | sed -n ${{median_line}}p > {wildcards.family}_{wildcards.sample}.median
        rm {wildcards.family}_{wildcards.sample}.dedup.cov
        """
