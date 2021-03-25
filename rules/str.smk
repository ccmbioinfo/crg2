rule EH:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EH/{sample}-{unit}.json",
        vcf = "str/EH/{sample}-{unit}.vcf",
        bam = "str/EH/{sample}-{unit}_realigned.bam"
    params:
        ref = config["ref"]["genome"],
        sex = lambda w: "`sh {}/scripts/str_helper.sh mapped/{}-{}.sorted.bam`".format(workflow.basedir, w.sample, w.unit),
        catalog = config["annotation"]["eh"]["catalog"]
    log:
        "logs/str/{sample}-{unit}-EH.log"
    wrapper:
        get_wrapper_path("EH")
