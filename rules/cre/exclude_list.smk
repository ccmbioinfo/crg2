rule validate_metadata:
    input:
        metadata=config["run"]["metadata"],
        samples=config["run"]["samples"]
    output:
        validated="excludelists/{family}.metadata.validated.tsv"
    log:
        "logs/metadata/{family}.validate_metadata.log"
    conda:
        "../../envs/crg.yaml"
    script:
        "../../scripts/validate_metadata.py"

rule write_sample_groups:
    input:
        metadata="excludelists/{family}.metadata.validated.tsv"
    output:
        baseline="excludelists/{family}.baseline_samples.txt",
        nonbaseline="excludelists/{family}.nonbaseline_samples.txt"
    log:
        "logs/metadata/{family}.write_sample_groups.log"
    conda:
        "../../envs/crg.yaml"
    script:
        "../../scripts/write_sample_groups.py"

rule baseline_germline_excludelist:
    input:
        vcf="annotated/gatk_haplotype/vep/{family}.gatk_haplotype.vep.vcf.gz",
        tbi="annotated/gatk_haplotype/vep/{family}.gatk_haplotype.vep.vcf.gz.tbi"
    output:
        vcf="excludelists/{family}.baseline_germline.vcf.gz"
    params:
        baseline_samples=get_baseline_samples(),
        min_dp=config["variant_filters"]["germline_exclude_build"]["min_dp"],
        min_alt_depth=config["variant_filters"]["germline_exclude_build"]["min_alt_depth"]
    log:
        "logs/metadata/{family}.baseline_germline_excludelist.log"
    conda:
        "../../envs/crg.yaml"
    script:
        "../../scripts/baseline_germline_excludelist.py"

rule exclude_baseline_germline_from_germline:
    input:
        vcf="annotated/gatk_haplotype/vep/{family}.gatk_haplotype.vep.vcf.gz",
        exclude="excludelists/{family}.baseline_germline.vcf.gz",
        vcf_tbi="annotated/gatk_haplotype/vep/{family}.gatk_haplotype.vep.vcf.gz.tbi",
        exclude_tbi="excludelists/{family}.baseline_germline.vcf.gz.tbi"
    output:
        vcf="study_filtered/gatk_haplotype/{family}.gatk_haplotype.baseline_excluded.vep.vcf.gz"
    log:
        "logs/metadata/{family}.exclude_baseline_germline_from_germline.log"
    conda:
        "../../envs/crg.yaml"
    shell:
        """
        bcftools isec -C -w1 -Oz -o {output.vcf} {input.vcf} {input.exclude} 2> {log}
        """

rule exclude_baseline_germline_from_somatic:
    input:
        vcf="annotated/gatk_mutect2/vep/{family}.gatk_mutect2.vep.vcf.gz",
        exclude="excludelists/{family}.baseline_germline.vcf.gz",
        vcf_tbi="annotated/gatk_mutect2/vep/{family}.gatk_mutect2.vep.vcf.gz.tbi",
        exclude_tbi="excludelists/{family}.baseline_germline.vcf.gz.tbi"
    output:
        vcf="study_filtered/gatk_mutect2/{family}.gatk_mutect2.baseline_excluded.vep.vcf.gz"
    log:
        "logs/metadata/{family}.exclude_baseline_germline_from_somatic.log"
    conda:
        "../../envs/crg.yaml"
    shell:
        """
        bcftools isec -C -w1 -Oz -o {output.vcf} {input.vcf} {input.exclude} 2> {log}
        """

rule final_filter_germline:
    input:
        vcf="study_filtered/gatk_haplotype/{family}.gatk_haplotype.baseline_excluded.vep.vcf.gz",
        tbi="study_filtered/gatk_haplotype/{family}.gatk_haplotype.baseline_excluded.vep.vcf.gz.tbi"
    output:
        vcf="study_filtered/gatk_haplotype/{family}.gatk_haplotype.final_filtered.vep.vcf.gz"
    params:
        nonbaseline_samples=get_non_baseline_samples(),
        min_dp=config["variant_filters"]["germline_final"]["min_dp"],
        min_alt_depth=config["variant_filters"]["germline_final"]["min_alt_depth"]
    log:
        "logs/metadata/{family}.final_filter_germline.log"
    conda:
        "../../envs/crg.yaml"
    script:
        "../../scripts/filter_annotated_germline.py"

rule final_filter_somatic:
    input:
        vcf="study_filtered/gatk_mutect2/{family}.gatk_mutect2.baseline_excluded.vep.vcf.gz",
        tbi="study_filtered/gatk_mutect2/{family}.gatk_mutect2.baseline_excluded.vep.vcf.gz.tbi"
    output:
        vcf="study_filtered/gatk_mutect2/{family}.gatk_mutect2.final_filtered.vep.vcf.gz"
    params:
        nonbaseline_samples=get_non_baseline_samples(),
        min_dp=config["variant_filters"]["somatic_final"]["min_dp"],
        min_af=config["variant_filters"]["somatic_final"]["min_af"],
        min_alt_ac=config["variant_filters"]["somatic_final"]["min_alt_ac"]
    log:
        "logs/metadata/{family}.final_filter_somatic.log"
    conda:
        "../../envs/crg.yaml"
    script:
        "../../scripts/filter_annotated_somatic.py"

rule per_sample_variant_counts:
    input:
        germline_annotated="annotated/gatk_haplotype/vep/{family}.gatk_haplotype.vep.vcf.gz",
        germline_baseline_excluded="study_filtered/gatk_haplotype/{family}.gatk_haplotype.baseline_excluded.vep.vcf.gz",
        germline_final="study_filtered/gatk_haplotype/{family}.gatk_haplotype.final_filtered.vep.vcf.gz",
        somatic_annotated="annotated/gatk_mutect2/vep/{family}.gatk_mutect2.vep.vcf.gz",
        somatic_baseline_excluded="study_filtered/gatk_mutect2/{family}.gatk_mutect2.baseline_excluded.vep.vcf.gz",
        somatic_final="study_filtered/gatk_mutect2/{family}.gatk_mutect2.final_filtered.vep.vcf.gz"
    output:
        summary="study_filtered/{family}.per_sample_variant_counts.tsv"
    params:
        germline_min_dp=config["variant_filters"]["germline_final"]["min_dp"],
        germline_min_alt_depth=config["variant_filters"]["germline_final"]["min_alt_depth"],
        somatic_min_dp=config["variant_filters"]["somatic_final"]["min_dp"],
        somatic_min_af=config["variant_filters"]["somatic_final"]["min_af"],
        somatic_min_alt_ac=config["variant_filters"]["somatic_final"]["min_alt_ac"]
    log:
        "logs/metadata/{family}.per_sample_variant_counts.log"
    conda:
        "../../envs/crg.yaml"
    script:
        "../../scripts/per_sample_variant_counts.py"
