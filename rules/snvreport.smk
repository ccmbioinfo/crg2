rule allsnvreport:
    input:
        db="annotated/gemini.db",
        vcf="annotated/all.vep.vcfanno.pass.vcf"
    output:
        directory("report/all")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/all-cre.log"
    threads: 2
    resources:
         mem_mb=40000
    params:
         cre=config["cre"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p report/all/{project}
         ln -s ../../{input.db} report/all/{project}/{project}-ensemble.db
         bgzip {input.vcf} -c -@ 2 > report/all/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s ./{project}-gatk-haplotype-annotated-decomposed.vcf.gz report/all/{project}/{project}-ensemble-annotated-decomposed.vcf.gz
         cd report/all
         {params.cre}/cre.sh {project} {params.ref}
         '''

if "panel" in config["run"] and config["run"]["panel"]: #present in config file and non-empty string
    rule panelsnvreport:
        input:
            db="annotated/gemini.db",
            vcf="annotated/all.vep.vcfanno.pass.vcf",
            panel=config["run"]["panel"]
        output:
            directory("report/panel")
        conda:
            "../envs/cre.yaml"
        log:
            "logs/report/panel-cre.log"
        threads: 2
        resources:
            mem_mb=40000
        params:
            cre=config["cre"],
            ref=config["ref"]["genome"]
        shell:
             '''
             mkdir -p report/panel/{project}
             ln -s ../../{input.db} report/panel/{project}/{project}-ensemble.db
             bedtools intersect --header -a {input.vcf} -b {input.panel} | bgzip -c -@ 2 > report/panel/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz
             ln -s ./{project}-gatk-haplotype-annotated-decomposed.vcf.gz report/panel/{project}/{project}-ensemble-annotated-decomposed.vcf.gz
             cd report/panel
             {params.cre}/cre.sh {project} {params.ref}
             '''

    rule panelflanksnvreport:
        input:
            db="annotated/gemini.db",
            vcf="annotated/all.vep.vcfanno.pass.vcf",
            panel=config["run"]["panel"],
        output:
            dir=directory("report/panel-flank"),
            panelflank="panel-{params.flank}.bed"
        conda:
            "../envs/cre.yaml"
        log:
            "logs/report/panel-flank-cre.log"
        threads: 2
        resources:
            mem_mb=40000
        params:
            cre=config["cre"],
            ref=config["ref"]["genome"],
            flank=config["run"]["flank"]
        shell:
             '''
             mkdir -p {output.dir}/{project}
             ln -s ../../{input.db} {output.dir}/{project}-ensemble.db
             cd {output.dir}
             cat {input.panel} | awk -F "\t" '{{print $1"\t"$2-{params.flank}"\t"$3+{params.flank}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output.panelflank}
             bedtools intersect --header -a {input.vcf} -b {output.panelflank} | bgzip -c -@ 2 > {project}-gatk-haplotype-annotated-deco$
             ln -s ./{project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}-ensemble-annotated-decomposed.vcf.gz
             {params.cre}/cre.sh {project} {params.ref}
             '''
