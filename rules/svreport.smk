rule svreport:
    input:
        get_annotated_sv_vcf()
    output:
        dir=directory("report/sv")
    log:
        "logs/report/sv.log"
    params:
        hgmd=config["annotation"]["svreport"]["hgmd"],
        hpo=config["run"]["hpo"],
        protein_coding_genes=config["annotation"]["svreport"]["protein_coding_genes"],
        exon_bed=config["annotation"]["svreport"]["exon_bed"],
        exac=config["annotation"]["svreport"]["exac"],
        omim=config["annotation"]["svreport"]["omim"],
        gnomad=config["annotation"]["svreport"]["gnomad"],
        biomart=config["annotation"]["svreport"]["biomart"],
        mssng_manta_counts=config["annotation"]["svreport"]["mssng_manta_counts"],
        mssng_lumpy_counts=config["annotation"]["svreport"]["mssng_lumpy_counts"],
    shell:
        """
        mkdir -p {output.dir}
        cd {output.dir}
        python3 ~/crg/crg.intersect_sv_vcfs.py -protein_coding_genes={params.protein_coding_genes} -exon_bed={params.exon_bed} \
        -hgmd={params.hgmd} -hpo={params.hpo} -exac={params.exac} -omim={params.omim} -biomart={params.biomart} -gnomad={params.gnomad} \
        -sv_counts {params.mssng_manta_counts} ${params.mssng_lumpy_counts} \
        -o={{"%s.wgs.sv.v%s.%s.tsv" % (project, PIPELINE_VERSION, date.today().strftime("%Y-%m-%d"))}} -i {{",".join(["../%s" % vcf for vcf in input[0]])}}
        """
