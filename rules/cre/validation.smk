rule benchmark:
    input: 
        vcf="annotated/coding/{family}-ensemble-decomposed.vcf.gz",
        tbi="annotated/coding/{family}-ensemble-decomposed.vcf.gz.tbi",
        benchmark = config["validation"]["benchmark"],
        eval_bed = "{family}-callable-highconf.bed"
    output: 
        directory("validation/{family}"),
    log: 
        "logs/benchmark/{family}-vcfeval.log"
    params:
        java_opts=config["params"]["rtg-tools"]["java_opts"],
        pipeline=config["run"]["pipeline"],
        sdf = config["params"]["rtg-tools"]["vcfeval"]["sdf"]
    resources:
        mem=30
    wrapper:
        get_wrapper_path("rtg-tools", "vcfeval")

rule evaluation_bed:
    input:
        callable_bed = "mapped/{family}-sort-callable.bed",
        benchmark = config["validation"]["benchmark"]
    output:
        "{family}-callable-highconf.bed"
    log:
        "logs/benchmark/{family}-bed-merge.log"
    params:
        pipeline=config["run"]["pipeline"],
        family=config["run"]["project"]
    resources:
        mem=30
    shell:
        '''
        truth_bed=`grep -v "family" {input.benchmark} | grep "{params.family}" | grep "{params.pipeline}" | awk '{{ print $5;}}'`;
        echo ${{truth_bed}};
        bedtools intersect -a ${{truth_bed}} -b {input.callable_bed} | sort -k1,1 -k2,2n > {output}
        '''
