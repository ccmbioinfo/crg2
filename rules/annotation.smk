rule bakta:
    input:
        "shovill/{sra_run}"
    params:
        outdir="annotated/bakta/",
        db=config["tools"]["bakta_db"]
    output:
        directory("annotated/{sra_run}_bakta")
    log:
        "logs/bakta/{sra_run}.log"
    wrapper:
        get_wrapper_path("bakta")