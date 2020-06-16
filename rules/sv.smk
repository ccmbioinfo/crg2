rule manta:
    input:
        bam = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam",
        bai = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        outdir = directory("sv/manta/{sample}-{unit}/"),
        include_chroms = "--region " + \
            " --region ".join([c for c in get_contigs().tolist() if is_nonalt(c)])
    threads:
        4
    output:
        "sv/manta/{sample}-{unit}/results/variants/diploidSV.vcf.gz"
    log:
        "logs/manta/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("manta")

rule wham:
    input:
        bam = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam",
        bai = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        include_chroms = ",".join(
            [c for c in get_contigs().tolist() if is_nonalt(c)])
    threads:
        4
    output:
        "sv/wham/{sample}-{unit}.wham.vcf"
    log:
        "logs/wham/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("wham")


rule samblaster:
    # smoove can actually extract splitters and discordants on it's own,
    # but keeping this in here in case it's needed for something else
    input:
        bam = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam"
    output:
        splitters = "recal/{sample}-{unit}.splitters.bam",
        discordants = "recal/{sample}-{unit}.discordants.bam"
    wrapper:
        get_wrapper_path("samblaster")

rule smoove:
    input:
        bam = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam",
        bai = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        name = "sv/smoove/{sample}-{unit}",
        exclude_chroms = ",".join([c for c in get_contigs().tolist() if is_nonalt(c) == False])
    threads:
        4
    output:
        "sv/smoove/{sample}-{unit}-smoove.genotyped.vcf.gz"
    log:
        "logs/smoove/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("smoove")

rule metasv:
    input:
        bam = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam",
        bai = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"],
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
    log:
        "logs/metasv/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("metasv")


rule snpeff:
    input:
        "sv/metasv/{sample}-{unit}/variants.vcf.gz"
    output:
        "sv/metasv/{sample}-{unit}/variants.snpeff.vcf"
    log:
        "logs/snpeff/{sample}-{unit}.log"
    params:
        memory_initial = "-Xms750m",
        memory_max = "-Xmx20g",
        reference = "GRCh37.75",
        data_dir = "/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/snpeff"
    wrapper:
        get_wrapper_path("snpeff")
