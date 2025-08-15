rule prefetch_sra:
    input:
        accession=config["run"]["accession"]
    params:
        outdir = "sra/",
        sratoolkit = config["tools"]["sratoolkit"]
    output:
        sra = "sra/{sra_run}/{sra_run}.sra"
    shell:
        '''
        {params.sratoolkit}/prefetch --output-directory {params.outdir} {sra_run}
        '''

rule fasterq_dump:
    input:
        sra = "sra/{sra_run}/{sra_run}.sra"
    params:
        outdir = "fastq/",
        sratoolkit = config["tools"]["sratoolkit"]
    output:
        fastq1 = "fastq/{sra_run}_1.fastq",
        fastq2 = "fastq/{sra_run}_2.fastq",
    shell:
        '''
        {params.sratoolkit}/fasterq-dump --outdir {params.outdir} {input.sra}
        '''
