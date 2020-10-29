rule vcf_to_tsv:
    input:
        "annotated/vcfanno/all.vep.vcfanno.vcf"
    output:
        report("tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt -g --fmt DP AD --info ANN < {input} | "
        "gzip > {output}"


rule plot_stats:
    input:
        "tables/calls.tsv.gz"
    output:
        depths=report("plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
