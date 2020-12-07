def get_cre_bams(ext="bam"):
    return ["recal/{}-{}.{}".format(s,units.loc[s].unit[0],ext) for s in samples.index]

def get_cre_vcfs():
    return ["called/{}-{}-vt.vcf".format(config["run"]["project"],i) for i in ["freebayes", "platypus", "samtools_call"] ]
    

rule freebayes:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        #regions="/path/to/region-file.bed"
    output:
        "called/{project}-freebayes.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/{project}.log"
    params:
        extra=config["params"]["freebayes"],         # optional parameters
        chunksize=100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 8
    wrapper:
        get_wrapper_path("freebayes")