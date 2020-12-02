callers = ["freebayes", "platypus", "samtools_call"]

def get_aligned_bams(wildcards,ext):
    return expand("mapped/{wildcards.sample}-{unit}.ext", sample=wildcards.sample, unit=units.loc[wildcards.sample].unit)

rule freebayes:
    input:
        bam=get_aligned_bams("bam"),
        bai=get_aligned_bams("bam.bai"),
        ref=config["ref"]["genome"],
        #regions="/path/to/region-file.bed"
    output:
        "called/{project}.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/{project}.log"
    params:
        extra="",         # optional parameters
        chunksize=100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 8
    wrapper:
        get_wrapper_path("freebayes")

rule platypus:
    input:
	# single or list of bam files
        bam=get_aligned_bams("bam")
        bai=get_aligned_bams("bam.bai")
        ref=config["ref"]["genome"],
        regions="regions.bed" #remove or empty quotes if not using regions
    output:
	    "called/{project}-platypus.vcf"
    threads: 8
    log:
        "logs/platypus/{project}.log"
    wrapper:
       get_wrapper_path("platypus")

rule samtools_call:
    input:
        bam=get_aligned_bams("bam")
        bai=get_aligned_bams("bam.bai")
        ref=config["ref"]["genome"],
        #regions="regions.bed" #remove or empty quotes if not using regions
    output:
	    "called/{project}-samtools_call.vcf"
    threads: 8
    log:
        "logs/samtools_call/{project}.log"
    wrapper:
       get_wrapper_path("samtools", "call")

rule vt:
    input:
        expand("called/{project}-{caller}.vcf", caller=callers)
    output:
        expand("called/{project}-{caller}-vt.vcf", caller=callers)
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{project}-{caller}.log"
    wrapper:
        get_wrapper_path("vt")


rule pass:
    input:
        vcf=expand("called/{project}-{caller}-vt.vcf", caller=callers)
    output:
	    expand("called/{project}-{caller}-vt-pass.vcf", caller=callers)
    params: "-f PASS"
    threads: 8
    resources:
        mem_mb = 4000
    log:
        "logs/pass/{project}-{caller}.log"
    wrapper:
       get_wrapper_path("bcftools", "view")


rule vep:
    input:
        expand("called/{project}-{caller}-vt-pass.vcf", caller=callers)
    output:
	    expand("annotated/vep/{project}-{caller}-vt-pass-vep.vcf", caller=callers)
    threads: 8
    log:
        "logs/vep/{project}-{caller}.log"
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        maxentscan=config["annotation"]["vep"]["maxentscan"],
        human_ancestor_fasta=config["annotation"]["vep"]["human_ancestor_fasta"],
        ref=config["ref"]["genome"],
    wrapper:
       get_wrapper_path("vep")


#merge vcfs:  merge vcfs from  all callers -> table (use criteria to retain calls b/w samples + caller priority)
#rule vcfanno: 
#rule vcf2db
#reporting