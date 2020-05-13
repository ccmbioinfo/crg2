rule snpeff:
    input:
        "filtered/all.vcf.gz",
    output:
        calls=report("annotated/all.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats="snpeff/all.csv",
	stats="snpeff/stats.csv",
    log:
        "logs/snpeff.log",
    params:
        reference=config["ref"]["name"],
        extra="-Xmx6g",
	data_dir="/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/snpeff",
    wrapper:
        get_wrapper_path("snpeff")
