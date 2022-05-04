

rule freebayes_mosaic:
    input:
        samples=get_cre_bams(),
        bai=get_cre_bams(ext="bam.bai"),
        ref=config["ref"]["genome"],
        regions="mapped/bed/{family}-sort-callable-{{contig}}.bed".format(family=project)
    output:
        temp("called/freebayes_mosaic/{contig}.vcf") 
    log:
        "logs/freebayes/mosaic/{contig}.log"
    params:
        extra=config["params"]["freebayes"]["mosaic_call"], 
        chunksize=100000
    threads: 1
    resources:
        mem=lambda wildcards, threads: threads * 4 if threads > 1 else 20
    wrapper:
        get_wrapper_path("freebayes")
    
rule merge_freebayes_mosaic:
    input:
        expand("called/freebayes_mosaic/{contig}.vcf",contig = get_canon_contigs())
    output:
        protected("genotyped/{family}-freebayes_mosaic.vcf")
    log:
        "logs/freebayes/mosaic/{family}-merge.log"
    shell:
        '''
        bcftools concat {input} | bcftools sort > {output}  
        '''

# calls will be bgzipped and tabixed with rule bgzip and tabix from cre/calling.smk
