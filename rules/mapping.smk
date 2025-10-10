
rule input_prep:
    input:
        units=config["run"]["units"]
    params:
        outdir = temp("fastq/"),
        sort_check = True,
        ref_cache = config["ref"]["ref_cache"],
        orad = config["tools"]["orad"],
        orad_ref = config["ref"]["orad_ref"]
    output:
        fastq1 = temp("fastq/{family}_{sample}_R1.fastq.gz"),
        fastq2 = temp("fastq/{family}_{sample}_R2.fastq.gz")
    log:
        "logs/input_prep/{family}_{sample}.log"
    threads:
        4
    wrapper:
        get_wrapper_path("samtools", "fastq")
    
        
rule map_reads:
    input:
        reads = ["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        temp("mapped/{family}_{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{family}_{sample}.log"
    params:
        index = config["ref"]["genome"],
        extra = get_read_group,
        verbosity = config["params"]["bwa"]["verbosity"],
        maxMem = config["params"]["bwa"]["maxMem"],
        markSplitReads = config["params"]["bwa"]["markSplitReads"],
        sort = "samtools",
        sort_order = "coordinate"
    threads: 16
    resources: 
        mem=lambda wildcards, threads: threads * 2
    wrapper:
        get_wrapper_path("bwa", "mem")


rule mark_duplicates:
    input:
        "mapped/{family}_{sample}.sorted.bam"
    output:
        bam = temp("dedup/{family}_{sample}.bam"),
        metrics = "qc/dedup/{family}_{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{family}_{sample}.log"
    params:
        markDuplicates = config["params"]["picard"]["MarkDuplicates"],
        java_opts = config["params"]["picard"]["java_opts"],
    wrapper:
        get_wrapper_path("picard", "markduplicates")


rule recalibrate_base_qualities:
    input:
        bam = "dedup/{family}_{sample}.bam",
        bai = "dedup/{family}_{sample}.bam.bai",
        ref = config["ref"]["genome"],
        known = config["ref"]["known_variants"],
        bed = "dedup/{family}_{sample}-callable.bed" 
    output:
        bam = temp("recal/{family}_{sample}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/bqsr/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "baserecalibrator")


rule mosdepth:
    input:
        bam = "dedup/{family}_{sample}.bam",
        bai = "dedup/{family}_{sample}.bam.bai",
    output:
        qbed = "dedup/{family}_{sample}.quantized.bed.gz",
        bed = "dedup/{family}_{sample}-callable.bed"
    params:
        prefix = f"dedup/{{family}}_{{sample}}",
        by = config["ref"]["canon_bed"],
        quantiles = "0:1:4:"
    log:
        "logs/mosdepth/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("mosdepth")
       

rule md5:
    input: 
        cram = "recal/{family}_{sample}.cram"
    output:
        md5 = "recal/{family}_{sample}.cram.md5"
    shell:
        """
        md5sum {input.cram} > {output.md5}
        """

rule remove_decoy:
    #redirect first samtools command to a temp output file "decoy.bam" 
    #otherwise it floods the log file with binary stream of decoy reads
    input:
        bam = "recal/{family}_{sample}.bam",
        canon = config["ref"]["canon_bed"],
    output:
        out_f = temp("decoy_rm/{family}_{sample}.no_decoy_reads.bam")
    log:
        "logs/remove_decoys/{family}_{sample}.log"
    threads: 8
    resources: mem_mb = 10000
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view --threads {threads} -h -L {input.canon} {input.bam} | egrep -v "hs37d5|NC_007605" | samtools view --threads {threads} - -bh > {output.out_f}
        """

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        get_wrapper_path("samtools", "index")


rule concatenate_callable_bed:
    input:
        expand("dedup/{family}_{sample}-callable.bed",family=project,sample=samples.index)
    output: 
        bed = "dedup/{family}-concat-sort.bed"
    params:
        tempdir = temp("dedup/temp/")
    log:
        "logs/bash/{family}-callable-concat.log"
    shell:
        '''
            export TMPDIR={params.tempdir}
            if [ ! -d {params.tempdir} ]; then mkdir -p {params.tempdir}; fi
            cat {input} | sort -k1,1 -k2,2n > {output.bed} 2>{log}
        '''
    
rule merge_bed:
    input: 
        "dedup/{family}-concat-sort.bed"
    output:
        "dedup/{family}-sort-callable.bed"
    log:
        "logs/bedtools/{family}-callable-merge.log"
    wrapper:
        get_wrapper_path("bedtools", "merge")

rule contigwise_bed:
    input:
        "dedup/{family}-sort-callable.bed"
    output:
        "dedup/bed/{family}-sort-callable-{contig}.bed"
    log:
        "logs/bash/{family}.{contig}.log"    
    shell:
        """
            awk '{{ if($1=="{wildcards.contig}") print $0; }}' {input} > {output}
            # if there are no callable regions on MT, bed file is empty and freebayes throws an error
            if [ ! -s dedup/bed/{wildcards.family}-sort-callable-MT.bed ]; then echo -e "MT\t1\t10"  >  dedup/bed/{wildcards.family}-sort-callable-MT.bed; fi 2>{log}
        """

rule bam_to_cram:
    input:
        bam = "recal/{family}_{sample}.bam",
        ref = config["ref"]["genome"]
    output:
        protected("recal/{family}_{sample}.cram")
    threads: 8 
    log:
        "logs/samtools/view/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("samtools", "view")
    

rule samtools_index_cram:
    input:
        "{prefix}.cram"
    output:
        "{prefix}.cram.crai"
    wrapper:
        get_wrapper_path("samtools", "index")
