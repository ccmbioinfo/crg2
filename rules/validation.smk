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

rule bcftools_isec:
    input: 
        new_vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz",
        new_tbi="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz.tbi",
        old_vcf = config["params"]["bcftools"]["isec"]["old_vcf"],
        old_tbi = config["params"]["bcftools"]["isec"]["old_vcf"] + ".tbi",

    output:
        directory("isec/{family}")
        #expand("isec/{family}/000{n}.vcf.gz", n=range(4))
    params:
        outdir = "isec/{family}"
    shell:
        '''
        bcftools isec {input.new_vcf} {input.old_vcf} -p {params.outdir} -O z 
        for i in {params.outdir}/000*.vcf.gz; do 
        zgrep -v "#" $i |wc -l >> {params.outdir}/summary.txt; 
        done
        '''

#rule vcf_compare:
#    input: expand("isec/000{n}.vcf.gz", n=range(4))
#    output: directory("vcf_compare/{family}".format(family=config["run"]["project"])
#    params:	
#        family = config["run"]["project"]
#    shell:
#        '''
#            mkdir -p vcf_compare/{params.family}
#        '''
