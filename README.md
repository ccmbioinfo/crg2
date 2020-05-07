# crg2
Research pipeline for exploring clinically relevant genomic variants

crg2 is a research pipeline aimed at discovering clinially relevant variants (SNVs, SVs) in whole genome sequencing data.
It aims to provide reproducible results, be computationally efficient, and transparent in it's workflow.

crg2 uses Snakemake and Anaconda to manage jobs and software dependencies.

## Set up instructions

1. Download and setup Anaconda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
2. Install Snakemake via Conda: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
3. Git clone this repo
4. (TBD) Download the following data files and set their paths in config.yaml
5. (TBD) Set up vcfanno and link up the data files

## Running the pipeline
1. Make a new directory, copy over samples.tsv, units.tsv and edit them to point to fastq files as necessary.
2. To run the entire pipeline as a single job edit the paths of dnaseq.pbs to point to the Snakefile and the 
conda installation directory. Then, ```qsub dnaseq.pbs``` in the working directory.
3. Reports can be found in the report/{PROJECT_ID} directory

## Pipeline details

### SNV
Paired end reads are mapped to GRCh37d5. Duplicate reads are maked but not removed. GATK4 best practices are applied:
base recalibration, gVCF's are generated per sample, combined and joint genotyping is performed. VCF's are decomposed,
normalized, and uniq'd. VCF's are annotated with custom annotations using vcfanno. Then, VCF's are annotated using VEP.

The VCF is transformed into a gemini db using vcf2db.py. A report is generated using the cre pipeline and its custom filters.

### SV
TBD

## Report columns

SNV: https://docs.google.com/document/d/1zL4QoINtkUd15a0AK4WzxXoTWp2MRcuQ9l_P9-xSlS4

SV: https://docs.google.com/document/d/1o870tr0rcshoae_VkG1ZOoWNSAmorCZlhHDpZuZogYE
