rule bakta:
    input:
        "shovill/{sra_run}"
    params:
        outdir="annotated/{sra_run}_bakta/",
        db=config["tools"]["bakta_db"],
        user_proteins=None
    output:
        directory("annotated/{sra_run}_bakta")
    log:
        "logs/bakta/{sra_run}.log"
    wrapper:
        get_wrapper_path("bakta")

rule bakta_PI_90:
    input:
        "shovill/{sra_run}"
    params:
        outdir="annotated/{sra_run}_bakta_PI_90/",
        db=config["tools"]["bakta_db"],
        user_proteins="/hpf/projects/CTrost/ajain/bacterial_genomics_pipeline/data/pseudomonas_annot_protein_list_curated/curated_pseudomonas_proteins-bakta_short_format.fasta",
        prefix="{sra_run}_PI_90"
    output:
        directory("annotated/{sra_run}_bakta_PI_90")
    log:
        "logs/bakta/{sra_run}_PI_90.log"
    wrapper:
        get_wrapper_path("bakta")

rule bakta_PI_30:
    input:
        "shovill/{sra_run}"
    params:
        outdir="annotated/{sra_run}_bakta_PI_30/",
        db=config["tools"]["bakta_db"],
        user_proteins="/hpf/projects/CTrost/ajain/bacterial_genomics_pipeline/data/pseudomonas_annot_protein_list_curated/curated_pseudomonas_proteins-bakta_long_format-PI_30.fasta",
        prefix="{sra_run}_PI_30"
    output:
        directory("annotated/{sra_run}_bakta_PI_30")
    log:
        "logs/bakta/{sra_run}_PI_30.log"
    wrapper:
        get_wrapper_path("bakta")

rule bakta_PI_30_SQC_50:
    input:
        "shovill/{sra_run}"
    params:
        outdir="annotated/{sra_run}_bakta_PI_30_SQC_50/",
        db=config["tools"]["bakta_db"],
        user_proteins="/hpf/projects/CTrost/ajain/bacterial_genomics_pipeline/data/pseudomonas_annot_protein_list_curated/curated_pseudomonas_proteins-bakta_long_format-PI_30_SubjQuery_50.fasta",
        prefix="{sra_run}_PI_30_SQC_50"
    output:
        directory("annotated/{sra_run}_bakta_PI_30_SQC_50")
    log:
        "logs/bakta/{sra_run}_PI_30_SQC_50.log"
    wrapper:
        get_wrapper_path("bakta")


rule defensefinder:
    input:
        expand("shovill/{sra_run}",sra_run=project)
    params:
        defensefinder_models="/hpf/projects/CTrost/ajain/bacterial_genomics_pipeline/tools/defensefinder_models",
        input_type="nucleotide"
    output:
        directory("annotated/defensefinder/")
    log:
       expand("logs/defensefinder/{sra_run}_defensefinder.log",sra_run=project)
    wrapper:
        get_wrapper_path("defensefinder")

rule defensefinder_prot:
    input:
        bakta_dir=expand("annotated/{sra_run}_bakta",sra_run=project),
    params:
        defensefinder_models="/hpf/projects/CTrost/ajain/bacterial_genomics_pipeline/tools/defensefinder_models",
        input_type="protein",
        acc=project
    output:
        directory("annotated/defensefinder_ProtInput/")
    log:
       expand("logs/defensefinder_prot/{sra_run}_defensefinder.log",sra_run=project)
    wrapper:
        get_wrapper_path("defensefinder")