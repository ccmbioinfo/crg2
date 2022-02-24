#gatk VCF is missing because it is symlinked to ensemble VCF
callers = [ gatk + "_haplotype", "samtools", "freebayes", "platypus", gatk + "_somatic" ] 

rule allsnvreport:
    input:
        db="annotated/coding/{family}-gemini.db",
        vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz",
        caller_vcfs = expand("filtered/{family}-{caller}.uniq.normalized.decomposed.pass.vcf.gz", family=project, caller=callers)
    output:
        directory("report/coding/{family}")
    conda:
        "../../envs/cre.yaml"
    log:
        "logs/report/{family}/cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         database_path=config["annotation"]["cre"]["database_path"],
         ref=config["ref"]["genome"],
         family = config["run"]["project"]
    shell:
         '''
         mkdir -p {output}
         cd {output}
         ln -s ../../../{input.db} {params.family}-ensemble.db
         for i in {input.caller_vcfs}; do
            caller=`basename ${{i}} .vcf.gz | cut -d "." -f1 | cut -d "-" -f2`;
            if [[ "${{caller}}" == "gatk3_haplotype" ]] || [[ "${{caller}}" == "gatk_haplotype" ]]; then
                caller="gatk-haplotype";
            fi;
            ln -s ../../../${{i}} {params.family}-${{caller}}-annotated-decomposed.vcf.gz;
         done
         ln -s ../../../{input.vcf} {params.family}-ensemble-annotated-decomposed.vcf.gz
         tabix {params.family}-ensemble-annotated-decomposed.vcf.gz
         cd ../
         cre={params.cre} reference={params.ref} database={params.database_path} {params.cre}/cre.sh {params.family} 
         cre={params.cre} reference={params.ref} database={params.database_path} type=wes.synonymous {params.cre}/cre.sh {params.family}
         unset type
         '''


checkpoint qc_check:
    input:
        dups = [expand("qc/dedup/{family}_{sample}.metrics.txt", sample=samples.index,family=project)],
        alignment = "qc/multiqc/multiqc_data/multiqc_general_stats.txt"
    output:
        "qc/qual_check.txt"
    run:
        qual_dict = {}
        for dups in input['dups']:
            dups = pd.read_csv(dups, sep='\t', comment='#')
            dups_perc = dups["PERCENT_DUPLICATION"][0]
            if dups_perc >= 0.2:
                qual_dict['high_dup_percentage'] = [True]
                print(qual_dict['high_dup_percentage'])
                break
            else:
                qual_dict['high_dup_percentage'] = [False]
        for alignment in pd.read_csv(input['alignment'], sep='\t')['QualiMap_mqc-generalstats-qualimap-percentage_aligned'].values:
            if float(alignment) < 90:
                qual_dict['low_alignment_percentage'] = [True]
                break
            else:
                qual_dict['low_alignment_percentage'] = [False]
        qual_df = pd.DataFrame.from_dict(qual_dict, orient='columns')
        qual_df.to_csv("qc/qual_check.txt", sep='\t', index=False)


def check_dup(wildcards):
    qual_check = pd.read_csv(checkpoints.qc_check.get(**wildcards).output[0], sep='\t')
    dups_perc = qual_check['high_dup_percentage'].values[0]
    alignment = qual_check['low_alignment_percentage'].values[0]
    if dups_perc == True and alignment == True:
        return expand("coverage/{family}_{sample}/",sample=samples.index,family=project) +  ["report/coding/{family}"] + ["qc/multiqc/multiqc.html"]
    elif dups_perc == True:
        return expand("coverage/{family}_{sample}/",sample=samples.index,family=project) +  ["report/coding/{family}"]
    else:
        return ["report/coding/{family}"]


rule minio:
    input: 
        # if duplication rate is >20% in any sample, generate coverage metrics and copy these to MinIO
        check_dup
    output:
        directory(config['run']['minio'] + "/{family}")
    log:
        "logs/minio/{family}.log"
    run:
        import os 
        import glob
        import shutil
        outdir=output[0]
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for f in input:
            if "report" in f:
                reports = glob.glob(os.path.join(f, '*wes*'))
                for r in reports:
                    print(f"copying {r} to MinIO")
                    shutil.copy(r, outdir)
            elif "qc" in f:
                print(f"copying {r} to MinIO")
                shutil.copy(f, outdir)
            else:
                coverage_sample = os.path.join(outdir, f)
                coverage_parent = os.path.join(outdir, "coverage")
                if not os.path.exists(coverage_parent):
                    os.mkdir(coverage_parent)
                print("copying coverage reports to MinIO")
                shutil.copytree(f, coverage_sample, copy_function = shutil.copy)

            


