def get_baseline_pon_vcfs(wildcards):
    return expand(
        "called/gatk_mutect2_pon/baseline/{family}_{sample}.pon.vcf.gz",
        family=wildcards.family,
        sample=BASELINE_SAMPLES,
    )

def get_baseline_pon_vcf_indices(wildcards):
    return expand(
        "called/gatk_mutect2_pon/baseline/{family}_{sample}.pon.vcf.gz.tbi",
        family=wildcards.family,
        sample=BASELINE_SAMPLES,
    )

rule baseline_mutect2_for_pon:
    input:
        map="recal/{family}_{sample}.bam",
        fasta=config["ref"]["genome"],
        mgp_germline=config["ref"]["known_variants"],
        intervals=config["ref"]["mutect_intervals"]
    output:
        vcf="called/gatk_mutect2_pon/baseline/{family}_{sample}.pon.vcf.gz",
        stats="called/gatk_mutect2_pon/baseline/{family}_{sample}.pon.vcf.gz.stats"
    threads: 4
    log:
        "logs/gatk/mutect2_pon/baseline/{family}_{sample}.log"
    conda:
        "../../envs/pon.yaml"
    shell:
        """
        gatk Mutect2 \
            -R {input.fasta} \
            -I {input.map} \
            --germline-resource {input.mgp_germline} \
            -max-mnp-distance 0 \
            -O {output.vcf} \
            > {log} 2>&1
        """

rule pon_genomicsdb_import:
    input:
        vcfs=get_baseline_pon_vcfs,
        indices=get_baseline_pon_vcf_indices,
        fasta=config["ref"]["genome"],
        intervals=config["ref"]["mutect_intervals"]
    output:
        db=directory("called/gatk_mutect2_pon/panel_of_normals/{family}.pon_db")
    params:
        variants=lambda wildcards, input: " ".join("-V {}".format(vcf) for vcf in input.vcfs)
    log:
        "logs/gatk/mutect2_pon/panel_of_normals/{family}.genomicsdbimport.log"
    conda:
        "../../envs/pon.yaml"
    shell:
        """
        gatk GenomicsDBImport \
            -R {input.fasta} \
            -L {input.intervals} \
            --genomicsdb-workspace-path {output.db} \
            {params.variants} \
            > {log} 2>&1
        """

rule create_panel_of_normals:
    input:
        fasta=config["ref"]["genome"],
        db="called/gatk_mutect2_pon/panel_of_normals/{family}.pon_db"
    output:
        vcf="called/gatk_mutect2_pon/panel_of_normals/{family}.baseline.panel_of_normals.vcf.gz"
    log:
        "logs/gatk/mutect2_pon/panel_of_normals/{family}.log"
    conda:
        "../../envs/pon.yaml"
    shell:
        """
        gatk CreateSomaticPanelOfNormals \
            -R {input.fasta} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            > {log} 2>&1
        """

rule gatk_mutect2:
    input:
        map = "recal/{family}_{sample}.bam",
        fasta=config["ref"]["genome"],
        mgp_germline=config["ref"]["known_variants"],
        pon=get_pon_vcf
    output:
        vcf="called/gatk_mutect2/{family}_{sample}.vcf",
        stats="called/gatk_mutect2/{family}_{sample}.vcf.stats"
    threads: 4
    log:
        "logs/gatk/mutect2/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "mutect")


rule pileup_summaries:
    input:
        cram="recal/{family}_{sample}.cram",
        cram_index="recal/{family}_{sample}.cram.crai",
        common_variants=config["ref"]["known_variants_common"],
        ref=config["ref"]["genome"]
    output:
        "called/gatk_mutect2/{family}_{sample}.pileups.table"
    log:
        "logs/gatk/getpileupsummaries/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "getpileupsummaries")


rule calculate_contamination:
    input:
        pileups="called/gatk_mutect2/{family}_{sample}.pileups.table"
    output:
        segments="called/gatk_mutect2/{family}_{sample}.segments.table",
        table="called/gatk_mutect2/{family}_{sample}.contamination.table"
    log:
        "logs/gatk/calculatecontamination/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "calculatecontamination")
    

rule filter_mutect_call:
    input: 
        vcf="called/gatk_mutect2/{family}_{sample}.vcf",
        fasta=config["ref"]["genome"],
        contamination="called/gatk_mutect2/{family}_{sample}.contamination.table",
        segments="called/gatk_mutect2/{family}_{sample}.segments.table"
    output:
        vcf="genotyped/{family}_{sample}.gatk_somatic.vcf"
    log:
        "logs/gatk/filtermutectcalls/{family}_{sample}_somatic.log"
    params:
        extra=config["params"]["gatk"]["FilterMutectCalls"]
        #java_opts=config["params"]["gatk"]["java_opts"]
    wrapper:
        get_wrapper_path("gatk", "filtermutectcalls")


rule pass_somatic:
    input:
        "genotyped/{family}_{sample}.gatk_somatic.vcf.gz", "genotyped/{family}_{sample}.gatk_somatic.vcf.gz.tbi"
    output:
        "genotyped/{family}_{sample}.gatk_somatic.pass_mutect2.vcf.gz"
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: 
        samples = get_sample_order,
        filter = "-f 'PASS,.' "
    wrapper:
        get_wrapper_path("bcftools", "view")


rule merge_mutect2_sample:
    input:
        vcf=get_gatk_somatic_vcf(),
        index=get_gatk_somatic_vcf(ext="vcf.gz.tbi")
    output:
        vcf_unsort=temp("filtered/{family}.gatk_mutect2.unsorted.vcf"),
        vcf="filtered/{family}.gatk_mutect2.vcf"
    log: 
        "logs/bcftools/merge/{family}.gatk_mutect2.log"
    params:
        extra=config["params"]["bcftools"]["merge"]
    wrapper:
        get_wrapper_path("bcftools", "merge")
