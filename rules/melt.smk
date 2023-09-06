
rule melt_preprocess:
    input:                        
        bam = "recal/{family}_{sample}.bam",
        bam_bai = "recal/{family}_{sample}.bam.bai"
    output: 
        disc = temp("recal/{family}_{sample}.bam.disc"),
        disc_bai = temp("recal/{family}_{sample}.bam.disc.bai"),
        fq = temp("recal/{family}_{sample}.bam.fq")
    log:
        "logs/melt_preprocess/{family}_{sample}.log"
    params:
        melt = config["tools"]["melt"],
        genome_ref = config["ref"]["genome"]
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} Preprocess \
        -bamfile {input.bam} \
        -h {params.genome_ref}
        '''


rule melt_indiv_analysis:
    input: 
        bam = [expand("recal/{family}_{sample}.bam".format(family=family, sample=s)) for s in samples.index],
        bam_bai = [expand("recal/{family}_{sample}.bam.bai".format(family=family, sample=s)) for s in samples.index],
        disc = [expand("recal/{family}_{sample}.bam.disc".format(family=family, sample=s)) for s in samples.index],
        disc_bai = [expand("recal/{family}_{sample}.bam.disc.bai".format(family=family, sample=s)) for s in samples.index],
        fq = [expand("recal/{family}_{sample}.bam.fq".format(family=family, sample=s)) for s in samples.index]
    output: 
        indiv_dir = temp(directory("melt/IndivAnalysis_out/"))
    log:
        [expand("logs/melt_indiv_analysis/{family}_{sample}.log".format(family=family, sample=s)) for s in samples.index]
    params:
        melt = config["tools"]["melt"],
        genome_ref = config["ref"]["genome"],
        element_ref = config["ref"]["melt_element_ref"]
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} IndivAnalysis \
        -b hs37d5/NC_007605 \
        -c 40 \
        -h {params.genome_ref} \
        -bamfile {input.bam}\
        -t {params.element_ref}\
        -w melt/IndivAnalysis_out/
        '''


rule melt_group_analysis:
    input: 
        indiv_dir = "melt/IndivAnalysis_out/"
    output: 
        pre_geno_ALU = temp("melt/GroupAnalysis_out/ALU.pre_geno.tsv"),
        pre_geno_LINE1 = temp("melt/GroupAnalysis_out/LINE1.pre_geno.tsv"),
        pre_geno_SVA = temp("melt/GroupAnalysis_out/SVA.pre_geno.tsv")
    params:
        melt = config["tools"]["melt"],
        genome_ref = config["ref"]["genome"],
        family = family,
        element_ref = config["ref"]["melt_element_ref"],
        genes = config["annotation"]["melt"]["genes"]
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} GroupAnalysis \
        -h {params.genome_ref} \
        -n {params.genes} \
        -t {params.element_ref}\
        -discoverydir {input.indiv_dir} \
        -w melt/GroupAnalysis_out/ 
        '''
## Note: cannot add family name as prefix to output name as MELT "Genotype" step looks for specific output file name without any prefixes


rule melt_genotype:
    input: 
        pre_geno_ALU = "melt/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "melt/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "melt/GroupAnalysis_out/SVA.pre_geno.tsv",
        bam = "recal/{family}_{sample}.bam"
    output: 
        genotype_tsv_ALU = temp("melt/Genotype_out/{family}_{sample}.ALU.tsv"),
        genotype_tsv_LINE1 = temp("melt/Genotype_out/{family}_{sample}.LINE1.tsv"),
        genotype_tsv_SVA = temp("melt/Genotype_out/{family}_{sample}.SVA.tsv")
    log:
        "logs/melt_genotype/{family}_{sample}.log"
    params:
        melt = config["tools"]["melt"],
        genome_ref = config["ref"]["genome"],
        element_ref = config["ref"]["melt_element_ref"]
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} Genotype \
        -bamfile {input.bam} \
        -h {params.genome_ref} \
        -t {params.element_ref}\
	    -p melt/GroupAnalysis_out/ \
        -w melt/Genotype_out/;
        '''


rule melt_make_VCF:
    input: 
        pre_geno_ALU = "melt/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "melt/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "melt/GroupAnalysis_out/SVA.pre_geno.tsv",
        genotype_tsv_ALU =  [expand("melt/Genotype_out/{family}_{sample}.ALU.tsv".format(family=family, sample=s)) for s in samples.index],
        genotype_tsv_LINE1 =  [expand("melt/Genotype_out/{family}_{sample}.LINE1.tsv".format(family=family, sample=s)) for s in samples.index],
        genotype_tsv_SVA =  [expand("melt/Genotype_out/{family}_{sample}.SVA.tsv".format(family=family, sample=s)) for s in samples.index]
    output: 
        VCF_ALU = temp("melt/MakeVCF_out/{family}_ALU.final_comp.vcf"),
        VCF_LINE1 = temp("melt/MakeVCF_out/{family}_LINE1.final_comp.vcf"),
        VCF_SVA = temp("melt/MakeVCF_out/{family}_SVA.final_comp.vcf")
    log:
        "logs/melt_make_VCF/{family}.log"
    params:
        melt = config["tools"]["melt"],
        genome_ref = config["ref"]["genome"],
        family = family,
        element_ref = config["ref"]["melt_element_ref"]
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} MakeVCF \
        -genotypingdir melt/Genotype_out/ \
        -p melt/GroupAnalysis_out/ \
        -h {params.genome_ref} \
        -t {params.element_ref}\
        -w melt/MakeVCF_out/ ; 

        for me in "ALU" "LINE1" "SVA"; do 
            mv "$me".final_comp.vcf melt/MakeVCF_out/{wildcards.family}_"$me".final_comp.vcf;
            mv melt/MakeVCF_out/"$me".hum.list melt/MakeVCF_out/{wildcards.family}_"$me".hum.list
        done;

        '''


rule melt_merge_VCFs:
    input: 
        VCF_ALU = "melt/MakeVCF_out/{family}_ALU.final_comp.vcf.gz",
        VCF_LINE1 = "melt/MakeVCF_out/{family}_LINE1.final_comp.vcf.gz",
        VCF_SVA = "melt/MakeVCF_out/{family}_SVA.final_comp.vcf.gz",

        VCF_ALU_index = "melt/MakeVCF_out/{family}_ALU.final_comp.vcf.gz.tbi",
        VCF_LINE1_index = "melt/MakeVCF_out/{family}_LINE1.final_comp.vcf.gz.tbi",
        VCF_SVA_index = "melt/MakeVCF_out/{family}_SVA.final_comp.vcf.gz.tbi"
    output:
        merged_VCF = "melt/MakeVCF_out/{family}_MELT.vcf.gz"
    log:
        "logs/melt_merge_VCFs/{family}.log"
    params:
        family = family
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        bcftools concat -a -Oz  -o melt/MakeVCF_out/{wildcards.family}_MELT.vcf.gz {input.VCF_ALU} {input.VCF_LINE1} {input.VCF_SVA};
        '''


rule melt_filter_VCF:
    input: 
        merged_VCF = "melt/MakeVCF_out/{family}_MELT.vcf.gz",
        index = "melt/MakeVCF_out/{family}_MELT.vcf.gz.tbi"
    output:
        filtered_VCF = "melt/MakeVCF_out/{family}_MELT_no_ac0.vcf.gz"
    log:
        "logs/melt_filter_VCF/{family}.log"
    params:
        family = family
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        bcftools filter -i 'FILTER!="ac0"' {input.merged_VCF} -Oz -o melt/MakeVCF_out/{wildcards.family}_MELT_no_ac0.vcf.gz
        '''


rule melt_snpeff:
    input: 
        filtered_VCF = "melt/MakeVCF_out/{family}_MELT_no_ac0.vcf.gz"
    output: 
        vcf = "melt/snpeff_out/{family}_MELT_snpeff.vcf",
        report = "melt/snpeff_out/{family}_MELT_snpeff_summary.html"
    log:
        "logs/melt_snpeff/{family}.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["dataDir"]
    wrapper:
        get_wrapper_path("snpeff")



rule melt_report:
    input: 
        vcf = expand("melt/snpeff_out/{family}_MELT_snpeff.vcf".format(family=family))
    output: 
        report_dir = directory("report/MELT/")
    log:
        expand("logs/melt_report/{family}.log".format(family=family))
    params:
        family = family,
        hpo = config["run"]["hpo"],
        protein_coding_genes = config["annotation"]["svreport"]["protein_coding_genes"],
        hgmd_db = config["annotation"]["svreport"]["hgmd"],
        exon_bed = config["annotation"]["svreport"]["exon_bed"],
        exac = config["annotation"]["svreport"]["exac"],
        omim = config["annotation"]["svreport"]["omim"],
        gnomad = config["annotation"]["svreport"]["gnomad_ins"], 
        biomart = config["annotation"]["svreport"]["biomart"],
        mssng_manta_counts = config["annotation"]["svreport"]["mssng_manta_counts"]

    conda: 
        "../envs/crg.yaml"

    script:
        "../scripts/MELT_report.py"
