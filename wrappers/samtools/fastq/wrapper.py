from snakemake.shell import shell
import os
import pandas as pd

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

family = snakemake.wildcards.family
sample = snakemake.wildcards.sample

name = f"{snakemake.wildcards.family}_{snakemake.wildcards.sample}"

units = pd.read_table(snakemake.input.units, dtype=str).set_index(["sample"], drop=False)
cram_file = units.loc[sample, 'cram']

dest = os.path.join(f"fastq/{family}_{sample}.cram")

copy_cmd = (" cp {cram_file} {dest}; ")

get_header_cmd = (" samtools view -H {dest} > fastq/{family}_{sample}_header.sam; ")

sed_cmd = (" sed -i 's+{snakemake.params.OldCramRef}+UR:{snakemake.params.NewCramRef}+g' fastq/{family}_{sample}_header.sam; ")

reheader_cmd = (" samtools reheader -i fastq/{family}_{sample}_header.sam {dest}; ")

fastq_cmd = (" samtools sort -n {dest} | samtools fastq - --reference {snakemake.params.NewCramRef} -1 fastq/{family}_{sample}_R1.fastq.gz -2 fastq/{family}_{sample}_R2.fastq.gz -s fastq/{family}_{sample}_singleton.fastq.gz; ")

rm_cmd = (" rm -f fastq/{family}_{sample}_header.sam fastq/{family}_{sample}.cram fastq/{family}_{sample}_singleton.fastq.gz; ")

shell("(" + copy_cmd + get_header_cmd + sed_cmd + reheader_cmd + fastq_cmd + rm_cmd + ") {log}")
