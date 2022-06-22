from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

java_opts = snakemake.params.get("java_opts")


shell(
    "export _JAVA_OPTIONS=\"{java_opts}\" && "
    "rtg vcfsubset -i {snakemake.input} -o annotated/{snakemake.wildcards.p}/vcfanno/noDP.vcf.gz --remove-info DP && "
    "bcftools +fill-tags annotated/{snakemake.wildcards.p}/vcfanno/noDP.vcf.gz -o {snakemake.output} -- -t 'DP2=sum(DP)' && "
    "sed -i \"s/DP2/DP/g\" {snakemake.output} && "
    "rm annotated/{snakemake.wildcards.p}/vcfanno/noDP.vcf* "

)