#callers = [ "freebayes_mosaic" ]
# {p} wildcard is ["freebayes_mosaic", "gatk_somatic"]
rule allsnvreport_mosaic:
    input:
        db="annotated/{p}/{family}-mosaic-gemini.db",
        vcf="annotated/{p}/vcfanno/{family}-mosaic.coding.vep.vcfanno.wDP.vcf.gz",
        caller_vcfs = "filtered/{family}-{p}.uniq.normalized.decomposed.vcf.gz"
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
        
        ln -s ../../../{input.caller_vcfs} {params.family}-mosaic-annotated-decomposed.vcf.gz
        
        ln -s ../../../{input.vcf} {params.family}-ensemble-annotated-decomposed.vcf.gz
        tabix {params.family}-ensemble-annotated-decomposed.vcf.gz
        cd ../
        {params.cre}/cre.sh {params.family} 
        type=wes.mosaic {params.cre}/cre.sh {params.family}
        unset type
        '''

if config["run"]["panel"]:

    def get_panel(wildcards):
                return "genes/{family}.bed"

    rule mosaic_panel2bed:
        input: 
            hpo=config["run"]["panel"], # hpo is equivalent to panel tsv file
            ensembl=config["genes"]["ensembl"],
            refseq=config["genes"]["refseq"],
            hgnc=config["genes"]["hgnc"]
        params: 
            crg2=config["tools"]["crg2"],
            cre=config["tools"]["cre"]
        output: 
            genes="genes/{family}.bed"
        conda: "../../envs/hpo_to_panel.yaml"
        log: 
            "logs/hpo_to_panel/{family}.log"
        script:
            "../../scripts/hpo_to_panel.py"



    rule mosaic_bamslice:
        input:
            genes=get_panel,
            bams=get_sample_bams
        output:
            bamslice = "genes/{family}_{sample}_slice.bam"
        params:
            extra = "panel"
        log: 
            "logs/mosaic_bamslice/{family}_{sample}.log"
        wrapper:
            get_wrapper_path("bedtools", "intersect")
            # need to modify the wrapper
            # bedtools intersect -abam ../recal/428_CH0034.bam -b 428.bed -ubam > 428_CH0034_test.bam


else:
    rule gatk_call_mosaic:
        input:
            map = "recal/{family}_{sample}.bam",
            #map="genes/{family}_{sample}_slice.bam",
            fasta=config["ref"]["genome"]
        output:
            vcf="called/gatk_mutect/{family}_{sample}.vcf",
            stats="called/gatk_mutect/{family}_{sample}.vcf.stats"
            #vcf=temp("called/gatk/{family}_{sample}.{contig}.g.vcf.gz")
        log:
            "logs/gatk/mutect/{family}_{sample}.log"
        params:
            extra=config["params"]["gatk"]["Mutect"], # in config this field is ""
            #java_opts=config["params"]["gatk"]["java_opts"]
        wrapper:
            get_wrapper_path("gatk", "mutect")


    rule filter_mutect_call:
        input: 
            vcf="called/gatk_mutect/{family}_{sample}.vcf",
            fasta=config["ref"]["genome"]
        output:
            vcf="genotyped/gatk_mutect/{family}_{sample}_somatic.vcf"
        log:
            "logs/gatk/filtermutectcalls/{family}_{sample}_somatic.log"
        params:
            extra=config["params"]["gatk"]["FilterMutectCalls"]
            #java_opts=config["params"]["gatk"]["java_opts"]
        wrapper:
            get_wrapper_path("gatk", "filtermutectcalls")
    
    rule merge_mutect_sample:
        input:
            vcf=get_gatk_somatic_vcf(),
            index=get_gatk_somatic_vcf(ext="vcf.gz.tbi")
        output:
            vcf_unsort=temp("genotyped/{family}-gatk_somatic_unsorted.vcf.gz"),
            vcf="genotyped/{family}-gatk_somatic.vcf"
        log: 
            "logs/bcftools/merge/{family}_gatk_somatic.log"
        params:
            extra=config["params"]["bcftools"]["merge"]
        wrapper:
            get_wrapper_path("bcftools", "merge")


    # filtering: skip softfilter step to do rule pass. Changed file name to distinguish from  "uniq.normalized.decomposed.pass.vcf.gz"
    rule pass_mosaic:
        input:
            "filtered/{prefix}.uniq.normalized.decomposed.vcf.gz", "filtered/{prefix}.uniq.normalized.decomposed.vcf.gz.tbi"
        output:
            "filtered/{prefix}_pass.uniq.normalized.decomposed.vcf.gz"
        threads: 6
        resources:
            mem=lambda wildcards, threads: threads * 2
        params: 
            samples = get_sample_order,
            filter = "-f 'PASS,.' "
        wrapper:
            get_wrapper_path("bcftools", "view")

    # make report
    # caller = "gatk_somatic"