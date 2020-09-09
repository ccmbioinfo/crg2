rule allsnvreport:
    input:
        db="annotated/gemini.db",
        vcf="annotated/vcfanno/all.vep.vcfanno.vcf"
    output:
        directory("report/all")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/all-cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p report/all/{project}
         ln -s ../../../{input.db} report/all/{project}/{project}-ensemble.db
         bgzip {input.vcf} -c > report/all/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s ./{project}-gatk-haplotype-annotated-decomposed.vcf.gz report/all/{project}/{project}-ensemble-annotated-decomposed.vcf.gz
         cd report/all
         {params.cre}/cre.sh {project}
         '''

if config["run"]["panel"]: #non-empty string
    rule panelsnvreport:
        input:
            db="annotated/gemini.db",
            vcf="annotated/vcfanno/all.vep.vcfanno.vcf",
            panel=config["run"]["panel"]
        output:
            dir=directory("report/panel")
        conda:
            "../envs/cre.yaml"
        log:
            "logs/report/panel-cre.log"
        resources:
            mem_mb=40000
        params:
            cre=config["tools"]["cre"],
            ref=config["ref"]["genome"],
        shell:
             '''
             mkdir -p {output.dir}/{project}
             ln -s ../../../{input.db} {output.dir}/{project}/{project}-ensemble.db
             bedtools intersect -header -a {input.vcf} -b {input.panel} | bgzip -c > {output.dir}/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz
             ln -s {output.dir}/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz {output.dir}/{project}/{project}-ensemble-annotated-decomposed.vcf.gz
             cd {output.dir}
             {params.cre}/cre.sh {project}
             '''

    rule panelflanksnvreport:
        input:
            db="annotated/gemini.db",
            vcf="annotated/vcfanno/all.vep.vcfanno.vcf",
            panel=config["run"]["panel"],
        output:
            dir=directory("report/panel-flank-{flank}"),
            panelflank="panel-{flank}.bed"
        conda:
            "../envs/cre.yaml"
        log:
            "logs/report/panel-flank-{flank}-cre.log"
        resources:
            mem_mb=40000
        params:
            cre=config["tools"]["cre"],
            ref=config["ref"]["genome"],
        shell:
             '''
             mkdir -p {output.dir}/{project}
             ln -s ../../../{input.db} {output.dir}/{project}/{project}-ensemble.db
             cat {input.panel} | awk -F "\t" '{{print $1"\t"$2-{flank}"\t"$3+{flank}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output.panelflank}
             bedtools intersect -header -a {input.vcf} -b {output.panelflank} | bgzip -c > {output.dir}/{project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz
             cd {output.dir}
             ln -s {project}/{project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}/{project}-ensemble-annotated-decomposed.vcf.gz
             {params.cre}/cre.sh {project}
             '''
