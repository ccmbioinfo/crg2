
callers = [ gatk + "_haplotype", "samtools", "freebayes", "platypus", gatk + "_somatic", "freebayes_mosaic" ] 
#callers = [ gatk + "_haplotype" ]

#list used for annotating VCFs with INFO/CALLERS
if len(callers) > 1:
    caller_annotation = [ "gatk-haplotype", "samtools", "freebayes", "platypus", "gatk-somatic", "freebayes-mosaic" ] 

def get_cre_vcfs():    
    return ["filtered/{family}-{caller}.uniq.normalized.decomposed.pass.vcf.gz".format(family=project,caller=i) for i in callers ]

def get_cre_vcf_tbi():
    return ["filtered/{family}-{caller}.uniq.normalized.decomposed.pass.vcf.gz.tbi".format(family=project,caller=i) for i in callers ]
    

rule vt:
    input:
        "genotyped/{prefix}.vcf.gz", "genotyped/{prefix}.vcf.gz.tbi"
    output:
        "filtered/{prefix}.uniq.normalized.decomposed.vcf"  
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{prefix}.log"
    wrapper:
        get_wrapper_path("vt")

rule soft_filter:
    input: 
        "filtered/{prefix}.uniq.normalized.decomposed.vcf.gz", "filtered/{prefix}.uniq.normalized.decomposed.vcf.gz.tbi" 
    output:
        "filtered/{prefix}.uniq.normalized.decomposed.softfiltered.vcf.gz" 
    params: 
        soft = config["filtering"]["soft"],
    log:
        "logs/bcftools/soft/{prefix}.log"
    wrapper:
        get_wrapper_path("bcftools", "filter")

rule pass:
    input:
        "filtered/{prefix}.uniq.normalized.decomposed.softfiltered.vcf.gz", "filtered/{prefix}.uniq.normalized.decomposed.softfiltered.vcf.gz.tbi"
    output:
        "filtered/{prefix}.uniq.normalized.decomposed.pass.vcf.gz"
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: 
        samples = get_sample_order,
        filter = "-f 'PASS,.' "
    wrapper:
        get_wrapper_path("bcftools", "view")

if len(get_cre_vcfs()) > 1:
    rule vcf_isec:
        input:
            vcf =  get_cre_vcfs(),
            tbi = get_cre_vcf_tbi()
        output:
            vcf = expand("isec/000{index}.vcf.gz", index=range(len(get_cre_vcfs()))),
            sites = "isec/sites.txt"
        params:
            outdir = "isec",
            filters = 'PASS',
            numpass = "+1",
        threads: 8
        log: "logs/isec.log"
        wrapper:
            get_wrapper_path("bcftools","isec")


    rule annot_caller:
        input: 
            all_sites = "isec/sites.txt",
            isec_vcf = expand("isec/000{index}.vcf.gz", index=range(len(get_cre_vcfs())))
        output: 
            all_sites = "isec/sites.caller.txt",
            callerwise_sites = expand("isec/{caller}.annot.{ext}", caller=caller_annotation,ext=["txt","txt.gz","txt.gz.tbi"]),
            callerwise_vcf = expand("isec/{caller}.annot.vcf.gz",caller=caller_annotation),
            hdr = "isec/hdr.txt"
        params:
            crg2 = config["tools"]["crg2"],
            callers = caller_annotation
        conda:
            "../../envs/common.yaml"
        shell:
            '''
            cat {input.all_sites} | parallel -k -j 16  {params.crg2}/scripts/annotate-caller.sh {{}} >> {output.all_sites}
            echo -e '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Variant called by"\\n##INFO=<ID=NUMCALLS,Number=1,Type=Integer,Description="Number of callers at this location">' > {output.hdr}
            for i in {params.callers}; do 
                sh {params.crg2}/scripts/callerwise_annotation.sh ${{i}} {output.all_sites} isec isec/${{i}}.annot.vcf.gz
            done;
            '''
    rule vcf_concat:
        input:
            vcf = expand("isec/{caller}.annot.vcf.gz",caller=caller_annotation),
            tbi = expand("isec/{caller}.annot.vcf.gz.tbi",caller=caller_annotation)

        output: 
            "filtered/{family}.numpass1.uniq.normalized.decomposed.vcf.gz"
        params: 
            "-a -d none"
        log: 
            "logs/bcftools/concat/{family}.log"
        threads: 8
        wrapper:
            get_wrapper_path("bcftools", "concat")

    rule ensemble:
        input: 
            "filtered/{family}.numpass1.uniq.normalized.decomposed.vcf.gz".format(family=project)
        output: 
            "annotated/coding/{family}.ensemble.decomposed.vcf.gz"
        log: 
            "logs/bcftools/filter/{family}.ensemble.log"
        params: 
            hard = " -i '(INFO/CALLERS=\"gatk-haplotype\" || INFO/NUMCALLS>=2)' "
        wrapper:
            get_wrapper_path("bcftools", "filter")
            

else:
    rule annot_caller:
        input: get_cre_vcfs()
        output: 
            txt = "isec/sites.caller.txt",
            bz = "isec/sites.caller.txt.gz",
            hdr = "isec/hdr.txt"
        conda:
            "../../envs/common.yaml"
        shell:
            '''
            mkdir -p isec
            bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\tgatk-haplotype\\n' {input} > {output.txt}
            bgzip -c {output.txt} > {output.bz}
            tabix -s1 -b2 -e2 {output.bz}
            echo -e '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Variant called by"\\n##INFO=<ID=NUMCALLS,Number=1,Type=Integer,Description="Number of callers at this location">' > {output.hdr}
            '''

    rule vcf_annotate:
        input: 
            vcf = get_cre_vcfs(),
            annot = "isec/sites.caller.txt.gz",
            hdr = "isec/hdr.txt"
        output: 
            "annotated/coding/{family}.ensemble.decomposed.vcf.gz"
        log: 
            "logs/bcftools/annotate/{family}.log"
        wrapper:
            get_wrapper_path("bcftools", "annotate")

