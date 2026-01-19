rule bakta:
    input:
        "shovill/{sra_run}"
    params:
        outdir="annotated/{sra_run}_bakta/",
        db=config["tools"]["bakta_db"]
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
        user_proteins="/hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/pseudomonas_annot_protein_list_curated/curated_pseudomonas_proteins-bakta_short_format.fasta",
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
        user_proteins="/hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/pseudomonas_annot_protein_list_curated/curated_pseudomonas_proteins-bakta_long_format-PI_30.fasta",
        prefix="{sra_run}_PI_30"
    output:
        directory("annotated/{sra_run}_bakta_PI_30")
    log:
        "logs/bakta/{sra_run}_PI_30.log"
    wrapper:
        get_wrapper_path("bakta")