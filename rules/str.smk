rule EH:
    input:
        bam="decoy_rm/{family}_{sample}.no_decoy_reads.bam",
        bai="decoy_rm/{family}_{sample}.no_decoy_reads.bam.bai"
    output:
        json = "str/EH/{family}_{sample}.json",
        vcf = "str/EH/{family}_{sample}.vcf",
        bam = "str/EH/{family}_{sample}_realigned.bam"
    params:
        ref = config["ref"]["genome"],
        sex = lambda w: "`sh {}/scripts/str/str_helper.sh decoy_rm/{}_{}.no_decoy_reads.bam`".format(workflow.basedir, w.family, w.sample),
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
        """
        echo "generating multi-sample genotypes" >> {log}
        python {params.crg2}/scripts/str/generate_EH_genotype_table.generic.py str/EH > {output.tsv}
        echo "annotating gene name & size threshold" >> {log}
        python {params.crg2}/scripts/str/add_gene+threshold_to_EH_column_headings2.py {output.tsv} {params.trf} > {output.annot}
        echo "generating final xlsx file" >> {log}
        python {params.crg2}/scripts/str/eh_sample_report.py {output.annot} {params.g1000} {output.xlsx} 
        prefix=`echo {output.xlsx} | awk '{{split($1,a,".xlsx"); print a[1]; }}'`;
        d=`date +%Y-%m-%d`
        outfile="${{prefix}}.${{d}}.xlsx";
        echo "Copying final report to filaname with timestamp: $outfile" >> {log}
        cp {output.xlsx} $outfile
        """

rule EHdn:
    input: expand("decoy_rm/{family}_{sample}.no_decoy_reads.bam",family=config["run"]["project"], sample=samples.index)
    output:
        json = "str/EHDN/{family}_EHDN_str.tsv"
    params:
        crg2 = config["tools"]["crg2"],
        ehdn = config["tools"]["ehdn"],
        g1k_manifest = config["annotation"]["ehdn"]["g1k_manifest"],
        ref = config["ref"]["genome"],
        # ref = config["ref"]["genome"],
        # irr_mapq = config["params"]["EHDN"]["irr_mapq"],
        # anchor_mapq = config["params"]["EHDN"]["anchor_mapq"],
    log:
        "logs/str/{family}-EHdn.log"
    conda:
        "../envs/ehdn.yaml"
    shell:
        '''
        EHDN={params.ehdn} g1k_manifest={params.g1k_manifest} ref={params.ref} script_dir={params.crg2}/scripts/str sh {params.crg2}/scripts/str/crg.ehdn.sh {wildcards.family} crg2
        '''

#rule EHdn_report:
#    input: "str/EHDN/{family}_EHDN_str.tsv".format(family=config["run"]["project"])
#    output: "report/str/{family}.EHDN.xlsx"
#    log: "logs/str/{family}_EHdn_report.log"
#    params:
#        repdir = "report/str",
#        crg2 = config["tools"]["crg2"],
#        family = config["run"]["project"],
#        outdir = "str/EHDN",
#    conda: "../envs/ehdn-report.yaml"
#    shell:
#        '''
#        export PATH="/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/bin/:$PATH"
#	    sh {params.crg2}/scripts/str/ehdn_report.sh {params.family} {params.outdir}
#        date=`date +%Y-%m-%d`;
#        f={params.outdir}/outliers/{params.family}.EHDN.${{date}}.xlsx;
#        if [ -f $f ]; then 
#	    if [ ! -d {params.repdir} ]; then mkdir -p {params.repdir}; fi;
#    	mv $f {params.repdir}
#	    cp {params.repdir}/{params.family}.EHDN.${{date}}.xlsx {params.repdir}/{params.family}.EHDN.xlsx
#        else exit; fi;
#        '''

rule EHDN_mark_outliers:
    input: "str/EHDN/{family}_EHDN_str.tsv".format(family=config["run"]["project"])
    output: "str/EHDN/{family}_outliers.txt".format(family=config["run"]["project"])
    params:
        crg2 = config["tools"]["crg2"],
        g1k_outlier = config["annotation"]["ehdn"]["g1k_outlier"]
    conda:
        "../envs/ehdn.yaml"
    shell:
        """
        python {params.crg2}/scripts/str/find_outliers.py {input} {output}
        cat {output} {params.g1k_outlier} > temp && mv temp {output}
        """

rule EHDN_DBSCAN_outlier:
    input: 
        profile = "str/EHDN/{family}_EHDN_str.tsv".format(family=config["run"]["project"]),
        outliers = "str/EHDN/{family}_outliers.txt".format(family=config["run"]["project"]),
    output: 
        clean="str/EHDN/outliers/clean.samples.txt",
        expansions="str/EHDN/outliers/EHdn.expansions.tsv",
        tr="str/EHDN/outliers/samples.with.TRcount.txt"

    params: 
        outdir="str/EHDN/outliers",
        crg2=config["tools"]["crg2"],
        trf=config["annotation"]["ehdn"]["trf"],
        g1k_samples=config["annotation"]["ehdn"]["g1k_samples"]
    conda: "../envs/ehdn-dbscan.yaml"
    shell:
        """
        Rscript {params.crg2}/scripts/str/DBSCAN.EHdn.parallel.R --infile {input.profile} --outpath {params.outdir} --outlierlist {input.outliers} --a1000g {params.g1k_samples} --exp {output.expansions}
        echo "Finished running DBSCAN.EHdn.parallel.R";
        echo "`ls {params.outdir}`";
        
        """

rule EHDN_merge_expansions:
    input: "str/EHDN/outliers/EHdn.expansions.tsv", "str/EHDN/{family}_EHDN_str.tsv".format(family=config["run"]["project"])
    output:  "str/EHDN/outliers/merged.rare.expansions.tsv","str/EHDN/outliers/merged.rare.expansions.forannotation.tsv",
            "str/EHDN/outliers/map.TRF.EHdn.0.66.tsv", "str/EHDN/outliers/merged.ehdn.tsv", 
    params: 
        outdir="str/EHDN/outliers",
        crg2=config["tools"]["crg2"],
        trf=config["annotation"]["ehdn"]["trf"],
        g1k_samples=config["annotation"]["ehdn"]["g1k_samples"]
    conda: "../envs/ehdn-dbscan.yaml"
    shell:
        """
        Rscript {params.crg2}/scripts/str/mergeExpansions.R --rscript {params.crg2}/scripts/str/ExpansionAnalysisFunctions.R --ehdn {input[1]}  --trf {params.trf} --outlier {input[0]} --outpath {params.outdir}
        """


rule EHDN_annovar:
    input: 
        merged_exp = "str/EHDN/outliers/merged.rare.expansions.tsv",
        ehdn_exp = "str/EHDN/outliers/EHdn.expansions.tsv"
    output: 
        merged_rare_exp = "str/EHDN/outliers/merged.rare.EHdn.expansion.tsv",
        omim_out = "str/EHDN/outliers/merged.rare.EHdn.expansion-OMIM.hg19_multianno.txt",
        gnomad_out = "str/EHDN/outliers/merged.rare.EHdn.expansion-gnoMAD.hg19_multianno.txt",
        xlsx="report/str/{family}.EHDN.xlsx".format(family=config["run"]["project"])
    params:
        crg2 = config["tools"]["crg2"],
        g1k_manifest = config["annotation"]["ehdn"]["g1k_manifest"],
        annovar = config["tools"]["annovar"],
        annovar_db = config["annotation"]["ehdn"]["annovar_db"],
        omim = config["annotation"]["ehdn"]["omim"],
        gnomad = config["annotation"]["ehdn"]["gnomad"],
        prefix = "str/EHDN/merged.rare.EHdn.expansion"
    conda: "../envs/eh-report.yaml"
    shell:
        """
            python  {params.crg2}/scripts/str/format_for_annovar.py {input.ehdn_exp} {input.merged_exp} {params.g1k_manifest}
            {params.annovar}/table_annovar.pl {output.merged_rare_exp} {params.annovar_db} -buildver hg19 -outfile {params.prefix}"-OMIM" -remove --onetranscript --otherinfo -protocol refGene -operation gx -nastring . -xreffile {params.omim}
            {params.annovar}/table_annovar.pl {output.merged_rare_exp} {params.annovar_db} -buildver hg19 -outfile {params.prefix}"-gnoMAD" -remove --onetranscript --otherinfo -protocol refGene -operation gx -nastring . -xreffile {params.gnomad}
            python {params.crg2}/scripts/str/format_from_annovar.py {output.gnomad_out} {output.omim_out} {output.xlsx}
        """