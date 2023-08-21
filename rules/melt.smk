
rule melt_preprocess:
    input:                        
        bam = "recal/{project}_{sample}.bam",
        bam_bai = "recal/{project}_{sample}.bam.bai"
    output: 
        disc = temp("recal/{project}_{sample}.bam.disc"),
        disc_bai = temp("recal/{project}_{sample}.bam.disc.bai"),
        fq = temp("recal/{project}_{sample}.bam.fq")
    log:
        "logs/melt_preprocess/{project}_{sample}.log"
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
        bam = [expand("recal/{project}_{sample}.bam".format(project=project, sample=s)) for s in samples.index],
        bam_bai = [expand("recal/{project}_{sample}.bam.bai".format(project=project, sample=s)) for s in samples.index],
        disc = [expand("recal/{project}_{sample}.bam.disc".format(project=project, sample=s)) for s in samples.index],
        disc_bai = [expand("recal/{project}_{sample}.bam.disc.bai".format(project=project, sample=s)) for s in samples.index],
        fq = [expand("recal/{project}_{sample}.bam.fq".format(project=project, sample=s)) for s in samples.index]
    output: 
        indiv_dir = temp(directory("melt/IndivAnalysis_out/"))
    log:
        [expand("logs/melt_indiv_analysis/{project}_{sample}.log".format(project=project, sample=s)) for s in samples.index]
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
        project = project,
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
## Note: cannot add project name as prefix to output name as MELT "Genotype" step looks for specific output file name without any prefixes


rule melt_genotype:
    input: 
        pre_geno_ALU = "melt/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "melt/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "melt/GroupAnalysis_out/SVA.pre_geno.tsv",
        bam = "recal/{project}_{sample}.bam"
    output: 
        genotype_tsv_ALU = temp("melt/Genotype_out/{project}_{sample}.ALU.tsv"),
        genotype_tsv_LINE1 = temp("melt/Genotype_out/{project}_{sample}.LINE1.tsv"),
        genotype_tsv_SVA = temp("melt/Genotype_out/{project}_{sample}.SVA.tsv")
    log:
        "logs/melt_genotype/{project}_{sample}.log"
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
        -w melt/Genotype_out/
        '''


rule melt_make_VCF:
    input: 
        pre_geno_ALU = "melt/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "melt/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "melt/GroupAnalysis_out/SVA.pre_geno.tsv",
        genotype_tsv_ALU =  [expand("melt/Genotype_out/{project}_{sample}.ALU.tsv".format(project=project, sample=s)) for s in samples.index],
        genotype_tsv_LINE1 =  [expand("melt/Genotype_out/{project}_{sample}.LINE1.tsv".format(project=project, sample=s)) for s in samples.index],
        genotype_tsv_SVA =  [expand("melt/Genotype_out/{project}_{sample}.SVA.tsv".format(project=project, sample=s)) for s in samples.index]
    output: 
        VCF_ALU = temp("melt/MakeVCF_out/{project}_ALU.final_comp.vcf"),
        VCF_LINE1 = temp("melt/MakeVCF_out/{project}_LINE1.final_comp.vcf"),
        VCF_SVA = temp("melt/MakeVCF_out/{project}_SVA.final_comp.vcf")
    log:
        "logs/melt_make_VCF/{project}.log"
    params:
        melt = config["tools"]["melt"],
        genome_ref = config["ref"]["genome"],
        project = project,
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
            mv "$me".final_comp.vcf melt/MakeVCF_out/{wildcards.project}_"$me".final_comp.vcf;
            mv melt/MakeVCF_out/"$me".hum.list melt/MakeVCF_out/{wildcards.project}_"$me".hum.list
        done;
        '''


rule melt_merge_VCFs:
    input: 
        VCF_ALU = "melt/MakeVCF_out/{project}_ALU.final_comp.vcf.gz",
        VCF_LINE1 = "melt/MakeVCF_out/{project}_LINE1.final_comp.vcf.gz",
        VCF_SVA = "melt/MakeVCF_out/{project}_SVA.final_comp.vcf.gz",

        VCF_ALU_index = "melt/MakeVCF_out/{project}_ALU.final_comp.vcf.gz.tbi",
        VCF_LINE1_index = "melt/MakeVCF_out/{project}_LINE1.final_comp.vcf.gz.tbi",
        VCF_SVA_index = "melt/MakeVCF_out/{project}_SVA.final_comp.vcf.gz.tbi"
    output:
        merged_VCF = "melt/MakeVCF_out/{project}_MELT.vcf.gz"
    log:
        "logs/melt_merge_VCFs/{project}.log"
    params:
        project = project
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        bcftools concat -a -Oz  -o melt/MakeVCF_out/{wildcards.project}_MELT.vcf.gz {input.VCF_ALU} {input.VCF_LINE1} {input.VCF_SVA};
        '''


rule melt_filter_VCF:
    input: 
        merged_VCF = "melt/MakeVCF_out/{project}_MELT.vcf.gz",
        index = "melt/MakeVCF_out/{project}_MELT.vcf.gz.tbi"
    output:
        filtered_VCF = "melt/MakeVCF_out/{project}_MELT_no_ac0.vcf.gz"
    log:
        "logs/melt_filter_VCF/{project}.log"
    params:
        project = project
    conda: 
        "../envs/melt.yaml"
    shell:
        '''
        bcftools filter -i 'FILTER!="ac0"' {input.merged_VCF} -Oz -o melt/MakeVCF_out/{wildcards.project}_MELT_no_ac0.vcf.gz
        '''


rule melt_snpeff:
    input: 
        filtered_VCF = "melt/MakeVCF_out/{project}_MELT_no_ac0.vcf.gz"
    output: 
        vcf = "melt/snpeff_out/{project}_MELT_snpeff.vcf",
        report = "melt/snpeff_out/{project}_MELT_snpeff_summary.html"
    log:
        "logs/melt_snpeff/{project}.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["dataDir"]
    wrapper:
        get_wrapper_path("snpeff")



rule melt_report:
    input: 
        vcf = expand("melt/snpeff_out/{project}_MELT_snpeff.vcf".format(project=project))
    output: 
        report_dir = directory("report/MELT/")
    log:
        expand("logs/melt_report/{project}.log".format(project=project))
    params:
        project = project,
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