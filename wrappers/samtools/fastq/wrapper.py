from snakemake.shell import shell
import pysam
import os
import pandas as pd
import fastq_prep

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    family = snakemake.wildcards.family
    sample = snakemake.wildcards.sample

    name = f"{family}_{sample}"

    units = pd.read_table(snakemake.input.units, dtype=str).set_index(
        ["sample"], drop=False
    )

    # determining the sample input type based on units.tsv
    file_types = units.loc[[sample], ~units.columns.isin(["sample", "platform"])]
    input_type = file_types.columns[-file_types.isnull().any()].tolist()
    print(input_type)
    if input_type == ["bam"]:
        print("File type is BAM")
        bam_file = units.loc[sample, "bam"]
        sort_cmd = ""
        bamtofastq_cmd = " samtools fastq {bam_file} -1 fastq/{family}_{sample}_R1.fastq.gz -2 fastq/{family}_{sample}_R2.fastq.gz -s /dev/null; "

        if snakemake.params.sort_check:
            print("Checking sort_order")
            bam = pysam.AlignmentFile(bam_file, "rb")
            header = bam.header
            sort_order = header["HD"]["SO"]
            bam.close()

            if sort_order != "queryname":
                print("Sorting by queryname")
                sort_cmd = "samtools sort -n -@ {snakemake.threads} -T {snakemake.params.outdir}/{snakemake.wildcards.sample} {bam_file} |"
                bamtofastq_cmd = "samtools fastq -1 fastq/{family}_{sample}_R1.fastq.gz -2 fastq/{family}_{sample}_R2.fastq.gz -s /dev/null; "

        shell("(" + sort_cmd + bamtofastq_cmd + ")")

    elif input_type == ["cram"]:
        print("File type is CRAM")
        cram_file = units.loc[sample, "cram"]
        dest = os.path.join(f"fastq/{family}_{sample}.cram")

        copy_cmd = " cp {cram_file} {dest} {log}; "
        if 'resarchivezone' in cram_file: # CRAM is stored in iRods archive on hpf
            copy_cmd = " module load irods_client; iget {cram_file} {dest} {log}; "

        fastq_cmd = (
            "export REF_CACHE={snakemake.params.ref_cache}; export REF_PATH={snakemake.params.ref_cache}; "
            " samtools sort -n -T {snakemake.params.outdir}/{snakemake.wildcards.sample} {dest} |  samtools fastq - "
            " -1 fastq/{family}_{sample}_R1.fastq.gz"
            " -2 fastq/{family}_{sample}_R2.fastq.gz"
            " -s fastq/{family}_{sample}_singleton.fastq.gz {log}; "
        )

        rm_cmd = " rm -f fastq/{family}_{sample}.cram fastq/{family}_{sample}_singleton.fastq.gz  {log}; "

        shell(
            "("
            + copy_cmd
            + fastq_cmd
            + rm_cmd
            + ")"
        )

    elif input_type == ["fq1", "fq2"]:
        print("File type is fastq")
        orad=snakemake.params.orad
        orad_ref=snakemake.params.orad_ref
        fastq_prep.main(units, sample, family, orad, orad_ref)
