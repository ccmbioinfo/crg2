rule benchmark:
    input: 
        vcf="annotated/coding/vcfanno/{family}-ensemble-decomposed.vep.vcfanno.vcf.gz",
        tbi="annotated/coding/vcfanno/{family}-ensemble-decomposed.vep.vcfanno.vcf.gz.tbi",
        benchmark = config["validation"]["benchmark"]
    output: 
        directory("validation/{family}"),
    log: 
        "logs/benchmark/{family}-vcfeval.log"
    params:
        java_opts=config["params"]["rtg-tools"]["java_opts"],
        pipeline="wes",
        sdf = config["params"]["rtg-tools"]["vcfeval"]["sdf"]
    resources:
        mem=30
    wrapper:
        get_wrapper_path("rtg-tools", "vcfeval")
