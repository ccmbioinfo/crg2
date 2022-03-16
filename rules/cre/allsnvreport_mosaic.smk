#callers = [ "freebayes_mosaic" ]

rule allsnvreport_mosaic:
    input:
        db="annotated/coding/{family}_mosaic-gemini.db",
        vcf="annotated/coding/vcfanno/{family}_mosaic.coding.vep.vcfanno.vcf.gz",
        #caller_vcfs is subject to change once soft filters have been decided. No filters applied yet!
        caller_vcfs = "filtered/{family}-freebayes_mosaic.uniq.normalized.decomposed.vcf.gz".format(family=project)
    output:
        directory("report/coding/mosaic/{family}")
    conda:
        "../../envs/cre.yaml"
    log:
        "logs/report/{family}/cre.log"
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
        ln -s ../../../../{input.db} {params.family}_mosaic-ensemble.db
        
        ln -s ../../../../{input.caller_vcfs} {params.family}-freebayes_mosaic-annotated-decomposed.vcf.gz
        
        ln -s ../../../../{input.vcf} {params.family}_mosaic-ensemble-annotated-decomposed.vcf.gz
        tabix {params.family}_mosaic-ensemble-annotated-decomposed.vcf.gz
        cd ../
        {params.cre}/cre.sh {params.family} 
        type=wes.mosaic {params.cre}/cre.sh {params.family}
        unset type
        '''