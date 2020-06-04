rule manta:
    input:
        bam = "recal/{sample}-{unit}.bam",
        fasta = config["ref"]["genome"]
    params:
        outdir = directory("sv/manta/{sample}-{unit}/")
    threads:
        4
    output:
        "sv/manta/{sample}-{unit}/results/variants/diploidSV.vcf.gz"
    wrapper:
        get_wrapper_path("manta")

rule wham:
    input:
        bam = "recal/{sample}-{unit}.bam",
        fasta = config["ref"]["genome"]
    params:
        # include MT as well?
        include_chroms = ",".join(
            [c for c in get_contigs().tolist() if is_autosomal_or_sex(c)])
    threads:
        4
    output:
        "sv/wham/{sample}-{unit}.wham.vcf"
    wrapper:
        get_wrapper_path("wham")


rule samblaster:
    # smoove can actually extract splitters and discordants on it's own,
    # but keeping this in here in case it's needed for something else
    input:
        bam = "recal/{sample}-{unit}.bam"
    output:
        splitters = "recal/{sample}-{unit}.splitters.bam",
        discordants = "recal/{sample}-{unit}.discordants.bam"
    wrapper:
        get_wrapper_path("samblaster")

rule smoove:
    input:
        bam = "recal/{sample}-{unit}.bam",
        fasta = config["ref"]["genome"]
    params:
        name = "sv/smoove/{sample}-{unit}"
    threads:
        4
    output:
        "sv/smoove/{sample}-{unit}-smoove.genotyped.vcf.gz"
    wrapper:
        get_wrapper_path("smoove")

rule metasv:
    input:
        bam = "recal/{sample}-{unit}.bam",
        fasta = config["ref"]["genome"],
        manta = "sv/manta/{sample}-{unit}/results/variants/diploidSV.vcf.gz",
        wham = "sv/wham/{sample}-{unit}.wham.vcf",
        lumpy = "sv/smoove/{sample}-{unit}-smoove.genotyped.vcf.gz"
    params:
        sample = "{sample}-{unit}",
        outdir = "sv/metasv/{sample}-{unit}"
    threads:
        4
    output:
        "sv/metasv/{sample}-{unit}/variants.vcf.gz"
    wrapper:
        get_wrapper_path("metasv")
