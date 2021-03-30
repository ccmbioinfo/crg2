rule manta:
    input:
        bam = "decoy_rm/{family}_{sample}.no_decoy_reads.bam",
        bai = "decoy_rm/{family}_{sample}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        outdir = directory("sv/manta/{family}_{sample}/"),
        include_chroms = "--region " + \
            " --region ".join([c for c in get_contigs().tolist() if is_nonalt(c)])
    threads:
        4
    output:
        "sv/manta/{family}_{sample}/results/variants/diploidSV.vcf.gz"
    log:
        "logs/manta/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("manta")

rule wham:
    input:
        bam = "decoy_rm/{family}_{sample}.no_decoy_reads.bam",
        bai = "decoy_rm/{family}_{sample}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        include_chroms = ",".join(
            [c for c in get_contigs().tolist() if is_nonalt(c)])
    threads:
        4
    output:
        "sv/wham/{family}_{sample}.wham.vcf"
    log:
        "logs/wham/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("wham")


rule samblaster:
    # smoove can actually extract splitters and discordants on it's own,
    # but keeping this in here in case it's needed for something else
    input:
        bam = "decoy_rm/{family}_{sample}.no_decoy_reads.bam"
    output:
        splitters = "recal/{family}_{sample}.splitters.bam",
        discordants = "recal/{family}_{sample}.discordants.bam"
    wrapper:
        get_wrapper_path("samblaster")

rule smoove:
    #smoove creates a lot intermediary bam files and .histo files;
    #added param for --outdir and use in the shell to delete the intermediary files with 'rm'
    #changed --name param to include only name
    input:
        bam = "decoy_rm/{family}_{sample}.no_decoy_reads.bam",
        bai = "decoy_rm/{family}_{sample}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"]
    params:
        outdir = "sv/smoove/{family}_{sample}/",
        name = "{family}_{sample}",
        exclude_chroms = ",".join([c for c in get_contigs().tolist() if is_nonalt(c) == False])
    threads:
       4
    output:
        "sv/smoove/{family}_{sample}/{family}_{sample}-smoove.genotyped.vcf.gz"
    log:
        "logs/smoove/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("smoove")

rule metasv:
    input:
        bam = "decoy_rm/{family}_{sample}.no_decoy_reads.bam",
        bai = "decoy_rm/{family}_{sample}.no_decoy_reads.bam.bai",
        fasta = config["ref"]["no_decoy"],
        manta = "sv/manta/{family}_{sample}/results/variants/diploidSV.vcf.gz",
        wham = "sv/wham/{family}_{sample}.wham.vcf",
        lumpy = "sv/smoove/{family}_{sample}/{family}_{sample}-smoove.genotyped.vcf.gz"
    params:
        sample = "{family}_{sample}",
        outdir = "sv/metasv/{family}_{sample}"
    threads:
        4
    output:
        "sv/metasv/{family}_{sample}/variants.vcf.gz"
    log:
        "logs/metasv/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("metasv")


rule snpeff:
    input:
        "sv/metasv/{family}_{sample}/variants.vcf.gz"
    output:
        vcf = "sv/metasv/{family}_{sample}/variants.snpeff.vcf",
        report = report("sv/metasv/{family}_{sample}/snpEff_summary.html",caption="../report/snpeff.rst",category="SnpEff")
    log:
        "logs/snpeff/{family}_{sample}.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["dataDir"]
    wrapper:
        get_wrapper_path("snpeff")

rule svscore:
    input:
        "sv/metasv/{family}_{sample}/variants.snpeff.vcf"
    output:
        "sv/metasv/{family}_{sample}/variants.snpeff.svscore.vcf"
    log:
        "logs/svscore/{family}_{sample}.log"
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
