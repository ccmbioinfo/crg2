rule snvreport:
    input:
        db="annotated/gemini.db",
        vcf="annotated/all.vep.vcfanno.pass.vcf"
    output:
        directory("report")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/cre.log"
    threads: 2
    resources:
         mem_mb=40000
    params:
         cre=config["cre"],
         ref=config["ref"]["genome"]
    shell:
         """
         mkdir -p report/{project}
         ln -s ../../{input.db} report/{project}/{project}-ensemble.db
         bgzip {input.vcf} -c -@ 2 > report/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s ./{project}-gatk-haplotype-annotated-decomposed.vcf.gz report/{project}/{project}-ensemble-annotated-decomposed.vcf.gz
         cd report
         {params.cre}/cre.sh {project} {params.ref}
         """
