rule allsnvreport:
    input:
        db="annotated/{p}/{family}-gemini.db",
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz"
    output:
        directory("report/{p}/{family}")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/{p}/{family}.cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         database_path=config["annotation"]["cre"]["database_path"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p {output}
         cd {output}
         ln -s ../../../{input.db} {project}-ensemble.db
         #bgzip ../../../{input.vcf} -c > {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s ../../../{input.vcf} {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}-ensemble-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {project}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         if [ {wildcards.p} == "coding" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} {params.cre}/cre.sh {project} 
         elif [ {wildcards.p} == "denovo" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} type=denovo {params.cre}/cre.sh {project} 
         else
         cre={params.cre} reference={params.ref} database={params.database_path} type=wgs {params.cre}/cre.sh {project}
         unset type
         fi;
         '''
if config["run"]["hpo"]:

    def get_panel(wildcards):
        if not config["run"]["panel"]:
            if wildcards.p == "panel":
                return "genes/{family}.bed"
            else:
                return "genes/{family}_{p}.bed"
        return config["run"]["panel"]

    # def get_bed(wildcards):
    #     if wildcards.p == "panel-flank":
    #         return "genes/{family}_{p}.bed"
    #     return get_panel()


    rule hpo_to_panel:
        input: 
            hpo=config["run"]["hpo"],
            ensembl=config["genes"]["ensembl"],
            refseq=config["genes"]["refseq"],
            hgnc=config["genes"]["hgnc"]
        params: 
            crg2=config["tools"]["crg2"],
            cre=config["tools"]["cre"]
        output: 
            genes="genes/{family}.bed"
        conda: "../envs/hpo_to_panel.yaml"
        log: "logs/hpo_to_panel/{family}.log"
        script:
            "../scripts/hpo_to_panel.py"

    rule add_flank:
        input: "genes/{family}.bed"
        output: "genes/{family}_{p}.bed"
        params: config["run"]["flank"]
        shell:
            '''
            cat {input} | awk -F "\t" '{{print $1"\t"$2-{params}"\t"$3+{params}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output}
            '''

    rule intersect:
        input: 
            left="filtered/{family}.vcf.gz",
            right=get_panel
        output:
            vcf="filtered/{p}/{family}.{p}.vcf.gz"
        params:
            extra="-header"
        log: "logs/report/bedtools-{family}-{p}.log"
        wrapper:
            get_wrapper_path("bedtools", "intersect")
            
        
    rule annotate_hpo:
        input:
            reports=expand("report/{p}/{family}",p=["coding", "panel", "panel-flank"], family=project),
            hpo=config["run"]["hpo"]
        output: 
            directory("report/hpo_annotated")
        conda: 
            "../envs/hpo_to_panel.yaml"
        params: 
            crg2 = config["tools"]["crg2"]
        log: 
            "logs/hpo_annotation.log"
        shell:
            '''
                if [ ! -d {output} ]; then mkdir -p {output}; fi;
                for i in {input.reports}; do 
                    echo "dir: ${{i}}" >> {log};                
                    if [[ "${{i}}" =~ .*"panel-flank".* ]]; then 
                        j=`find ${{i}} -name  "*.wgs.[0-9]*.csv" | grep -v "clinical"`;
                        rename=`echo ${{j}} | sed 's/wgs/wgs.panel-flank100k/g'`;
                        if [ ! -f ${{rename}} ]; then
                            ln -s `basename ${{j}}` ${{rename}};
                        else
                            echo "${{j}} found panel-flank"  >> {log};
                        fi;
                    elif [[ "${{i}}" =~ .*"panel".* ]]; then 
                        j=`find ${{i}} -name "*.wgs.[0-9]*.csv" | grep -v "clinical"`;
                        rename=`echo ${{j}} | sed 's/wgs/wgs.panel/g'`;
                        if [ ! -f ${{rename}} ]; then
                            ln -s `basename ${{j}}` ${{rename}};
                        else
                            echo "${{j}} found panel"  >> {log};
                        fi;
                    else 
                        rename=`find ${{i}} -name "*.wes*.[0-9]*.csv" | grep -v "clinical"`;
                        echo "${{rename}} found wes"  >> {log};
                    fi;
                    python {params.crg2}/scripts/add_hpo_terms_to_wes.py {input.hpo} ${{rename}} {output} >> {log} 2>&1
                done;
            '''