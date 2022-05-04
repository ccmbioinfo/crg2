rule vep_mosaic:
    input:
        "filtered/{family}-freebayes_mosaic.uniq.normalized.decomposed.vcf.gz"
    output:
        temp("annotated/mosaic/vep/{family}_mosaic.coding.vep.vcf")
    log:
        "logs/vep/{family}_mosaic.vep.coding.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        ref=config["ref"]["genome"],
    wrapper:
        get_wrapper_path("vep")

rule vcfanno_mosaic:
    input:
        "annotated/mosaic/vep/{family}_mosaic.coding.vep.vcf"
    output:
        "annotated/mosaic/vcfanno/{family}_mosaic.coding.vep.vcfanno.vcf"
    log:
        "logs/vcfanno/{family}_mosaic.vcfanno.coding.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
        lua_script = config["annotation"]["cre.vcfanno"]["lua_script"],
       	conf = config["annotation"]["cre.vcfanno"]["conf"],
        base_path = config["annotation"]["cre.vcfanno"]["base_path"],
    wrapper:
        get_wrapper_path("vcfanno")

rule fix_dp_mosaic:
    input:
        "annotated/mosaic/vcfanno/{family}_mosaic.coding.vep.vcfanno.vcf"
    output:
        "annotated/mosaic/vcfanno/{family}_mosaic.coding.vep.vcfanno.wDP.vcf"
    threads: 1
    resources:
        mem_mb = 20000
    shell:
        '''
        export _JAVA_OPTIONS=\"-Xmx2048m\"
        rtg vcfsubset -i {input} -o annotated/mosaic/vcfanno/noDP.vcf.gz --remove-info DP
        module load bcftools/1.11
        bcftools +fill-tags annotated/mosaic/vcfanno/noDP.vcf.gz -o {output} -- -t 'DP2=sum(DP)'
        sed -i "s/DP2/DP/g" {output}
        rm annotated/mosaic/vcfanno/noDP.vcf*
        '''

rule vcf2db_mosaic:
    input:
        "annotated/mosaic/vcfanno/{family}_mosaic.coding.vep.vcfanno.wDP.vcf".format(family=project)
    output:
         db = "annotated/mosaic/{family}_mosaic-gemini.db"
    log:
        "logs/vcf2db/vcf2db.{family}_mosaic.log"
    params:
        ped = config["run"]["ped"],
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")