#callers = [ "freebayes_mosaic" ]

rule allsnvreport_mosaic:
    input:
        db="annotated/mosaic/{family}_mosaic-gemini.db",
        vcf="annotated/mosaic/vcfanno/{family}_mosaic.coding.vep.vcfanno.wDP.vcf.gz",
        #caller_vcfs is subject to change once soft filters have been decided. No filters applied yet!
        caller_vcfs = "filtered/{family}-freebayes_mosaic.uniq.normalized.decomposed.vcf.gz".format(family=project)
    output:
        directory("report/mosaic/{family}")
    conda:
        "../../envs/cre.yaml"
    log:
        "logs/report/{family}/cre_mosaic.log"
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
        
        ln -s ../../../{input.caller_vcfs} {params.family}-freebayes_mosaic-annotated-decomposed.vcf.gz
        
        ln -s ../../../{input.vcf} {params.family}-ensemble-annotated-decomposed.vcf.gz
        tabix {params.family}-ensemble-annotated-decomposed.vcf.gz
        cd ../
        {params.cre}/cre.sh {params.family} 
        type=wes.mosaic {params.cre}/cre.sh {params.family}
        unset type
        '''