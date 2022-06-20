#gatk VCF is missing because it is symlinked to ensemble VCF
callers = [ gatk + "_haplotype", "samtools", "freebayes", "platypus" ] 

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

checkpoint dup_perc:
    input:
        [expand("qc/dedup/{family}_{sample}.metrics.txt", sample=samples.index,family=project)]
    output:
        "qc/dup_check.txt"
    run:
        for dups in input:
            dups = pd.read_csv(dups, sep='\t', comment='#')
            dups_perc = dups["PERCENT_DUPLICATION"][0]
            if dups_perc >= 0.2:
                with open(output[0], 'w') as f:
                    f.write('TRUE')
                    break
            else:
                with open(output[0], 'w') as f:
                    f.write('FALSE')


def check_dup(wildcards):
    with checkpoints.dup_perc.get(**wildcards).output[0].open() as f:
        if f.read().strip() == "TRUE":
            return expand("coverage/{family}_{sample}/",sample=samples.index,family=project) +  ["report/coding/{family}"]
        else:
            return ["report/coding/{family}"]


rule minio:
    input: 
        # if duplication rate is >20% in any sample, generate coverage metrics and copy these to MinIO
        check_dup
    output:
        directory("minio/{family}")
    log:
        "logs/minio/{family}.log"
    run:
        import os 
        import glob
        import shutil
        outdir=output[0]
        os.mkdir(outdir)
        for f in input:
            if "report" in f:
                reports = glob.glob(os.path.join(f, '*wes*'))
                for r in reports:
                    print(f"copying {r} to MinIO")
                    shutil.copy(r, outdir)
            else:
                print("copying coverage reports to MinIO")
                shutil.copytree(f, os.path.join(outdir, "coverage"), copy_function = shutil.copy)

            


