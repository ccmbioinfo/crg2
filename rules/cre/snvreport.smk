rule allsnvreport:
    input:
        db="annotated/{project}-gemini.db",
        vcf="annotated/vcfanno/{project}-ensemble-decomposed.vep.vcfanno.vcf"
    output:
        directory("report/{project}")
    conda:
        "../../envs/cre.yaml"
    log:
        "logs/report/{project}/cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p {output}
         cd {output}
         ln -s ../../{input.db} {project}-ensemble.db
         bgzip ../../{input.vcf} -c > {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}-ensemble-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {project}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         {params.cre}/cre.sh {project} 
         type=wes.synonymous {params.cre}/cre.sh {project}
         unset type
         '''
            
        
