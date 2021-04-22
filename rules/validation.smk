rule benchmark:
    input: 
        vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz",
        tbi="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz.tbi",
    output: 
        directory("validation/{family}"),
    log: 
        "logs/benchmark/{family}-vcfeval.log"
    params:
        java_opts=config["params"]["rtg-tools"]["java_opts"],
        sdf=config["params"]["rtg-tools"]["vcfeval"]["sdf"],
        baseline=config["params"]["rtg-tools"]["vcfeval"]["baseline"],
        eval_bed=config["params"]["rtg-tools"]["vcfeval"]["wgs_bed"]
    resources:
        mem=30
    wrapper:
        get_wrapper_path("rtg-tools", "vcfeval")

rule bgzip:
    input: "{prefix}.vcf"
    output: "{prefix}.vcf.gz"
    conda: "../envs/validation.yaml"
    shell:
        '''
        bgzip -c {input} > {output}
        '''
rule tabix:
    input: "{prefix}.vcf.gz"
    output: "{prefix}.vcf.gz.tbi"
    conda: "../envs/validation.yaml"
    shell:
        '''
        tabix {input}
        '''