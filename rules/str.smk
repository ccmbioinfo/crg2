rule EH:
    input:
        bam="mapped/{family}_{sample}.sorted.bam",
        bai="mapped/{family}_{sample}.sorted.bam.bai"
    output:
        json = "str/EH/{family}_{sample}.json",
        vcf = "str/EH/{family}_{sample}.vcf",
        bam = "str/EH/{family}_{sample}_realigned.bam"
    params:
        ref = config["ref"]["genome"],
        sex = lambda w: "`sh {}/scripts/str_helper.sh mapped/{}_{}.sorted.bam`".format(workflow.basedir, w.family, w.sample),
        catalog = config["annotation"]["eh"]["catalog"]
    log:
        "logs/str/{family}_{sample}-EH.log"
    wrapper:
        get_wrapper_path("EH")

rule EH_report:
    input:
        json = get_eh_json
    output:
        tsv = "str/EH/{family}_EH_str.tsv",
        annot = "str/EH/{family}_EH_str.annot.tsv",
        xlsx = "report/str/{family}.EH-v1.1.xlsx"
    log:
        "logs/str/{family}-eh-report.log"
    params:
        trf = config["annotation"]["eh"]["trf"],
        crg2 = config["tools"]["crg2"],
        g1000 = config["annotation"]["eh"]["1000g"]
    conda:
        "../envs/eh-report.yaml"
    shell:
        '''
        echo "generating multi-sample genotypes" > {log}
        python {params.crg2}/scripts/generate_EH_genotype_table.generic.py str/EH > {output.tsv}
        echo "annotating gene name & size threshold" > {log}
        python {params.crg2}/scripts/add_gene+threshold_to_EH_column_headings2.py {output.tsv} {params.trf} > {output.annot}
        echo "generating final xlsx file" > {log}
        python {params.crg2}/scripts/eh_sample_report.py {output.annot} {params.g1000} {output.xlsx} 
        prefix=`echo {output.xlsx} | awk '{{split($1,a,".xlsx"); print a[1]; }}'`;
        d=`date +%Y-%m-%d`
        outfile="${{prefix}}.${{d}}.xlsx";
        echo "Copying final report to filaname with timestamp: $outfile" > {log}
        cp {output.xlsx} $outfile
        '''
rule EHdn:
    input: expand("mapped/{family}_{sample}.sorted.bam",family=config["run"]["project"], sample=samples.index)
    output:
        json = "str/EHDN/{family}_EHDN_str.tsv"
    params:
        crg = config["tools"]["crg"],
        # ref = config["ref"]["genome"],
        # irr_mapq = config["params"]["EHDN"]["irr_mapq"],
        # anchor_mapq = config["params"]["EHDN"]["anchor_mapq"],
    log:
        "logs/str/{family}-EHdn.log"
    shell:
        '''
        sh {params.crg}/crg.ehdn.sh {wildcards.family} crg2
        '''
rule EHdn_report:
    input: "str/EHDN/{family}_EHDN_str.tsv".format(family=config["run"]["project"])
    output: directory("report/str")
    log: "logs/str/EHdn_report.log"
    params:
        crg = config["tools"]["crg"],
        family = config["run"]["project"],
        outdir = "str/EHDN"
    shell:
        '''
        sh {params.crg}/ehdn_report.sh {params.family} {params.outdir}
        date=`date +%Y-%m-%d`;
        f={params.outdir}/outliers/{params.family}.EHDN.${{date}}.xlsx;
        if [ -f $f ]; then mkdir -p report/str; mv $f report/str/ 
        else exit; fi;
        '''