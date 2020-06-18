# crg2
Clinical research pipeline for exploring variants in whole genome sequencing data

<div align="center">
    <img src="/crg2logolarge.png" width="800px"</img> 
</div>

crg2 is a research pipeline aimed at discovering clinically relevant variants (SNVs, SVs) in whole genome sequencing data.
It aims to provide reproducible results, be computationally efficient, and transparent in it's workflow.

crg2 uses Snakemake and Conda to manage jobs and software dependencies.

## Installation instructions

1. Download and setup Anaconda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
2. Install Snakemake 5.10.0 via Conda: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
3. Git clone this repo
4. Make a directory for Conda to install all it's environments and executables in, for example:
```
mkdir ~/crg2-conda
```

5. Navigate to the crg2 directory. Install all software dependencies using:
```
cd crg2
snakemake --use-conda -s Snakefile --conda-prefix ~/crg2-conda --create-envs-only
```
Make sure to replace ```~/crg2-conda``` with the path made in step 4. This will take a while.

6. Install these plugins for VEP: ```LoF, MaxEntScan, SpliceRegion```. Refer to this page for installation instructions: https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html. The INSTALL.pl script has been renamed to vep_install in the VEP's Conda build. It is located in the conda environment directory, under ```share/ensembl-vep-99.2-0/vep_install```. Therefore, your command should be similar to: ```fb5f2eb3/share/ensembl-vep-99.2-0/vep_install -a p --PLUGINS LoF,MaxEntScan,SpliceRegion```

7. Git clone cre: ```git clone https://github.com/ccmbioinfo/cre``` to a safe place

8. Replace the VEP path's to the VEP directory installed from step 6. Replace the cre path in crg2/config.yaml with the one from step 7.

## Running the pipeline
1. Make a folder in a directory with sufficient space. Copy over the template files samples.tsv, units.tsv, config.yaml.
You may need to re-copy config.yaml if the file was recently updated in repo for previous crg2 runs.
```
mkdir NA12878
cp crg2/samples.tsv crg2/units.tsv crg2/config.yaml NA12878
```

2. Reconfigure the 3 files to reflect project names, sample names, input fastq files, a panel bed file (if any) and a ped file (if any). Inclusion of a panel bed file will generate 2 SNV reports with all variants falling within these regions. Inclusion of a ped file currently does nothing except create a gemini db with the pedigree data stored in it.

samples.tsv
```
sample
NA12878
```

units.tsv
```
sample	unit	platform	fq1	fq2
NA12878	1	ILLUMINA	/hpf/largeprojects/ccm_dccforge/dccdipg/Common/NA12878/NA12878.bam_1.fq	/hpf/largeprojects/ccm_dccforge/dccdipg/Common/NA12878/NA12878.bam_2.fq
```

config.yaml
```
run:
  project: "NA12878"
  samples: samples.tsv
  units: units.tsv
  ped: "" # leave this line blank if there is no ped
  panel: "/hpf/largeprojects/ccmbio/dennis.kao/NA12878/panel.bed" # remove this line entirely if there is no panel bed file
  flank: "100000"
  
cre: /hpf/largeprojects/ccm_dccforge/dccdipg/Common/pipelines/cre

ref:
  name: GRCh37.75
  
...
```

3. Active the conda environment with Snakemake 5.10.0

```
(base) [dennis.kao@qlogin5 crg2]$ conda activate snakemake
(snakemake) [dennis.kao@qlogin5 crg2]$ snakemake -v
5.10.0
```

4. Test that the pipeline will run by adding the flag "-n" to the command in dnaseq.pbs. 

```
(snakemake) [dennis.kao@qlogin5 crg2]$ snakemake --use-conda -s /hpf/largeprojects/ccm_dccforge/dccdipg/Common/pipelines/crg2/Snakefile --cores 4 --conda-prefix ~/crg2-conda -n
Building DAG of jobs...
Job counts:
	count	jobs
	1	all
	86	call_variants
	86	combine_calls
	86	genotype_variants
	2	hard_filter_calls
	1	map_reads
	1	merge_calls
	1	merge_variants
	1	recalibrate_base_qualities
	2	select_calls
	1	snvreport
	1	vcf2db
	1	vcfanno
	1	vep
	1	vt
	272

[Tue May 12 10:56:12 2020]
rule map_reads:
    input: /hpf/largeprojects/ccm_dccforge/dccdipg/Common/NA12878/NA12878.bam_1.fq, /hpf/largeprojects/ccm_dccforge/dccdipg/Common/NA12878/NA12878.bam_2.fq
    output: mapped/NA12878-1.sorted.bam
    log: logs/bwa_mem/NA12878-1.log
    jobid: 271
    wildcards: sample=NA12878, unit=1
    threads: 4

[Tue May 12 10:56:12 2020]
rule recalibrate_base_qualities:
    input: mapped/NA12878-1.sorted.bam, /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa, /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/gemini_data/dbsnp.b147.20160601.tidy.vcf.gz
    output: recal/NA12878-1.bam
    log: logs/gatk/bqsr/NA12878-1.log
    jobid: 270
    wildcards: sample=NA12878, unit=1

[Tue May 12 10:56:12 2020]
rule call_variants:
    input: recal/NA12878-1.bam, /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa, /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/gemini_data/dbsnp.b147.20160601.tidy.vcf.gz
    output: called/NA12878.GL000204.1.g.vcf.gz
    log: logs/gatk/haplotypecaller/NA12878.GL000204.1.log
    jobid: 239
    wildcards: sample=NA12878, contig=GL000204.1
    
...
```

5. Run the pipeline as a single job using: ```qsub dnaseq.pbs```. The SNV report can be found in the directory: report/{PROJECT_ID}/{PROJECT_ID}.*.csv

## Pipeline details

### Pre-calling steps
1. Map fastq's to the human decoy genome GrCh37d5

2. Picard MarkDuplicates, but don't remove reads

3. GATK4 base recalibration

4. Remove reads mapped to decoy chromosomes

### SNV
1. Call SNV's and generate gVCFs

2. Merge gVCF's and perform joint genotyping

3. Filter against GATK best practices filters

4. Decompose multiallelics, sort and uniq the filtered VCF using vt

5. Annotate using vcfanno and VEP

6. Generate a gemini db using vcf2db.py

7. Generate a cre report using cre.sh

### SV
1. Call SV's using Manta, Smoove and Wham

2. Merge calls using MetaSV

3. Annotate VCF using snpEff and SVScores

4. Split multi-sample VCF into individual sample VCFs

5. Generate an annotated report using crg

## Reports

Column descriptions and more info on how variants are filtered can be found here:

SNV: https://docs.google.com/document/d/1zL4QoINtkUd15a0AK4WzxXoTWp2MRcuQ9l_P9-xSlS4

SV: https://docs.google.com/document/d/1o870tr0rcshoae_VkG1ZOoWNSAmorCZlhHDpZuZogYE

The pipeline generates 4 reports:

1. wgs.snv - a report on coding SNVs across the entire genome

2. wgs.panel.snv - a report on SNVs within the panel specified bed file

3. wgs.panel.snv - a report on SNVs within the panel specified bed file with a 100kb flank on each side

4. wgs.sv - a report on SVs across the entire genome
