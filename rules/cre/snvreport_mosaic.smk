# {p} wildcard is "gatk_somatic"
rule allsnvreport_mosaic:
    input:
        db="annotated/{p}/{family}-mosaic-gemini.db",
        vcf="annotated/{p}/vcfanno/{family}-mosaic.coding.vep.vcfanno.wDP.vcf.gz"
    output:
        directory("report/{p}/{family}")
    conda:
        "../../envs/cre.yaml"
    log:
        "logs/report/{p}/{family}.cre_mosaic.log"
    resources:
        mem_mb=40000
    params:
        cre=config["tools"]["cre"],
        ref=config["ref"]["genome"],
        family = config["run"]["project"]
    shell:
        '''
        mkdir -p {output}
        cd {output}
        ln -s ../../../{input.db} {params.family}-ensemble.db       
        ln -s ../../../{input.vcf} {params.family}-ensemble-annotated-decomposed.vcf.gz
        tabix {params.family}-ensemble-annotated-decomposed.vcf.gz
        cd ../
        cre={params.cre} type=wes.mosaic {params.cre}/cre.sh {params.family}
        '''
        