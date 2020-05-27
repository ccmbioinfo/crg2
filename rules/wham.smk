rule wham:
    input:
        bam = "/hpf/largeprojects/ccmbio/mcouse/testing/crg2/NA12878/recal/NA12878-1.bam",
        fasta = config["ref"]["genome"]
    params:
        # include MT as well?
        include_chroms = ",".join(
            [c for c in get_contigs().tolist() if is_autosomal_or_sex(c)])
    threads:
        4
    output:
        "sv/wham/wham.vcf"
    wrapper:
        get_wrapper_path("wham")
