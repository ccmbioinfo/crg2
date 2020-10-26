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
    #smoove creates a lot intermediary bam files and .histo files;
    #added param for --outdir and use in the shell to delete the intermediary files with 'rm'
    #changed --name param to include only name
    input:
        bam = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam",
        bai = "decoy_rm/{sample}-{unit}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        outdir = "sv/smoove",
        name = "{sample}-{unit}",
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
        temp("sv/metasv/{sample}-{unit}/variants.vcf.gz")
    log:
        "logs/metasv/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("metasv")


rule snpeff:
    input:
        "sv/metasv/{sample}-{unit}/variants.pass.vcf.gz"
    output:
        vcf = temp("sv/metasv/{sample}-{unit}/variants.snpeff.vcf"),
        report = report("sv/metasv/{sample}-{unit}/snpEff_summary.html",caption="../report/snpeff.rst",category="SnpEff")
    log:
        "logs/snpeff/{sample}-{unit}.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["dataDir"]
    wrapper:
        get_wrapper_path("snpeff")

rule svscore:
    input:
        "sv/metasv/{sample}-{unit}/variants.snpeff.vcf"
    output:
        "sv/metasv/{sample}-{unit}/variants.snpeff.svscore.vcf"
    log:
        "logs/svscore/{sample}-{unit}.log"
    params:
        svscore_script=config["tools"]["svscore_script"],
        operations=config["params"]["svscore"]["operations"],
        exon_bed=config["annotation"]["svscore"]["exon_bed"],
	intron_bed=config["annotation"]["svscore"]["intron_bed"],
	cadd=config["annotation"]["svscore"]["cadd"],
    conda:
        "../envs/svscore.yaml"
    shell:
        """
        perl -w {params.svscore_script} \
        -o {params.operations} \
        -e {params.exon_bed} \
        -f {params.intron_bed} \
        -c {params.cadd} \
        -i {input[0]} \
        > {output[0]}
        """
