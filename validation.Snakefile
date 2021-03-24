project = "NA12878"
wes = True
if wes:
    bed = "/hpf/largeprojects/ccmbio/ccmmarvin_shared/validation/benchmarking/results/NA12878-capture-callable-high_conf.bed"
else:
    bed = "/hpf/largeprojects/ccmbio/ccmmarvin_shared/validation/benchmarking/NISTv3.3.2/HG001_truth.bed"

rule benchmarking:
    input: 
        test = expand("annotated/{project}-ensemble-decomposed.vcf.{ext}", project=project, ext=("gz", "gz.tbi")),
        truth = "/hpf/largeprojects/ccmbio/ccmmarvin_shared/validation/benchmarking/NISTv3.3.2/HG001_truth.vcf.gz",
        bed = bed
    params:
        sdf = "/hpf/largeprojects/ccmbio/ccmmarvin_shared/validation/benchmarking/GRch37_SDF"
    output:
        directory("benchmarking")
    log:
        "logs/benchmarking.log"
    conda:
        "envs/rtg.yaml"
    shell:
        '''
        rtg vcfeval -b {input.truth} -c {input.test[0]} -t {params.sdf} -e {input.bed} -o {output}
        '''

        

