rule allsnvreport:
    input:
        db="annotated/{p}/gemini.db",
        vcf="annotated/{p}/vcfanno/all.{p}.vep.vcfanno.vcf"
    output:
        directory("report/{p}")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/{p}/cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p {output}/{project}
         cd {output}/{project}
         ln -s ../../../{input.db} {project}-ensemble.db
         bgzip ../../../{input.vcf} -c > {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}-ensemble-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {project}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         if [ {wildcards.p} == "coding" ]; then  
         {params.cre}/cre.sh {project} 
         else
         type=wgs {params.cre}/cre.sh {project}
         unset type
         fi;
         '''
if config["run"]["panel"]:

    def get_bed(wildcards):
        if wildcards.p == "panel-flank":
            return "{}-flank-{}k.bed".format(config["run"]["project"], int(config["run"]["flank"]/1000))
        else:
            return config["run"]["panel"]

    rule add_flank:
        input: config["run"]["panel"]
        output: "{}-flank-{}k.bed".format(config["run"]["project"], int(config["run"]["flank"]/1000))
        params: config["run"]["flank"]
        shell:
            '''
            cat {input} | awk -F "\t" '{{print $1"\t"$2-{params}"\t"$3+{params}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output}
            '''

    rule intersect:
        input: 
            left="filtered/all.vcf.gz",
            right=get_bed
        output:
            vcf="filtered/{p}/all.{p}.vcf.gz"
        params:
            extra="-header"
        log: "logs/report/bedtools-{p}.log"
        wrapper:
            get_wrapper_path("bedtools", "intersect")
            
        