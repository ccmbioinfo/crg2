rule EH:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EH/{sample}-{unit}.json",
        vcf = "str/EH/{sample}-{unit}.vcf",
        bam = "str/EH/{sample}-{unit}_realigned.bam"
    params:
        ref = config["ref"]["genome"],
        sex = lambda w: "`sh {}/scripts/str_helper.sh mapped/{}-{}.sorted.bam`".format(workflow.basedir, w.sample, w.unit),
        catalog = config["annotation"]["eh"]["catalog"]
    log:
        "logs/str/{sample}-{unit}-EH.log"
    wrapper:
        get_wrapper_path("EH")

rule EH_report:
    input:
        json = get_eh_json()
    output:
        tsv = "str/EH/{project}_EH_str.tsv",
        annot = "str/EH/{project}_EH_str.annot.tsv",
        xlsx = "report/str/{project}_EH_v1.1.xlsx"
    log:
        "logs/str/{project}-eh-report.log"
    params:
        trf = config["annotation"]["eh"]["trf"],
        crg2 = config["tools"]["crg2"],
        g1000 = config["annotation"]["eh"]["1000g"]
    conda:
        "../envs/eh-report.yaml"
    shell:
        '''
        python {params.crg2}/scripts/generate_EH_genotype_table.generic.py str/EH > {output.tsv} > {log} 2>&1
        python {params.crg2}/scripts/add_gene+threshold_to_EH_column_headings2.py {output.tsv} {params.trf} > {output.annot} > {log} 2>&1
        python {params.crg2}/scripts/eh_sample_report.py {output.annot} {params.g1000} {output.xlsx} > {log} 2>&1
        '''