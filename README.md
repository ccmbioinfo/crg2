# crg2-bacterial
Research pipeline for de-novo assembly of whole genome (WGS) of bacterial isolates. It aims to provide reproducible results, be computationally efficient, and transparent in it's workflow.crg2 uses Snakemake and Conda to manage jobs and software dependencies.

<div align="center">
    <img src="/crg2logolarge.png" width="800px"</img> 
</div>

## Conda environment set-up
Snakemake will build the conda-environment for each rule according to the environment.yaml file present in the respective wrapper's directory.

To set up the krona environment, I had to follow a few additional steps. Documenting it here for future reference:
1. Actiavte conda environment: `conda activate /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/0017e84a`
1. Remove the taxonomy directory:  `rm -rf /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/0017e84a/opt/krona/taxonomy`
2. Create the taxonomy directory in another folder:  `mkdir /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/conda_envs/krona/taxonomy`
3. Create a symlink: `ln -s /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/conda_envs/krona/taxonomy /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/0017e84a/opt/krona/taxonomy`
4. Change directory to the newly created taxonomy directory: `cd /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/conda_envs/krona/taxonomy`
5. Download the taxonomy file: `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
6. Build the taxonomy database:  `ktUpdateTaxonomy.sh --only-build`



## Pipeline details

### Pre-calling steps
1. Map fastqs to the human decoy genome GRCh37d5

2. Picard MarkDuplicates, but don't remove reads

3. GATK4 base recalibration

4. Remove reads mapped to decoy chromosomes

### WGS: SNV
1. Call SNV's and generate gVCFs

2. Merge gVCF's and perform joint genotyping

3. Filter against GATK best practices filters

4. Decompose multiallelics, sort and uniq the filtered VCF using vt

5. Annotate using vcfanno and VEP

6. Generate a gemini db using vcf2db.py

7. Generate a cre report using cre.sh

### WES: SNV
1. Call variants using GATK, Freebayes, Platypus, and SAMTools. 

2. Apply caller specific filters and retain PASS variants

3. Decompose multiallelics, sort and uniq filtered VCF using vt

4. Retain variants called by GATK or 2 other callers; Annotate caller info in VCF with INFO/CALLER and INFO/NUMCALLS.

5. Annotate using vcfanno and VEP

6. Generate a gemini db using vcf2db.py

7. Generate a cre report using cre.sh

8. Repeat the above 1-7 steps using GATK MUTECT2 variant calls to generate a mosaic variant report

### WGS: SV (filtered and unfiltered)
1. Call SV's using Manta, Smoove and Wham

2. Merge calls using MetaSV

3. Annotate VCF using snpEff and SVScores

4. Split multi-sample VCF into individual sample VCFs

5. Generate an annotated report (produces filtered, i.e. SV called by at least two callers, and unfiltered reports)

### WGS: BND
1. Filter Manta SV calls to exclude all variants but BNDs

2. Annotate VCF using snpEff

3. Generate an annotated BND report 


### WGS: STR

A. ExpansionHunter: known repeat location

  1. Identify repeat expansions in sample BAM/CRAMs
  2. Annotate repeats with disease threshold, gene name, repeat sizes from 1000Genome (mean& median) 
  3. Generate per-family report as Excel file

B. ExpansionHunterDenovo: denovo repeats

  1. Identify denovo repeat in sample BAM/CRAMs
  2. Combine individual JSONs from current family and 1000Genomes to a multi-sample TSV
  3. Run DBSCAN clustering to identify outlier repeats
  4. Annotate with gnoMAD, OMIM, ANNOVAR
  5. Generate per-family report as Excel file

### WGS: Mitochondrial variants

  1. Call MT variants using mity call
  2. Normalize MT variants using mity normalize
  3. Annotate MT variants using vcfanno
  4. Generate mitochondrial variant report

### WGS: MELT mobile element insertions (MEIs)
  1. Preprocess CRAMs
  2. Call MEIs in individuals
  3. Group and genotype MEIs across a family
  4. Filter out ac0 calls
  5. Annotate variants
  6. Generate MELT report

## Running the pipeline
1. Make a folder in a directory with sufficient space. Copy over the template files crg2/samples.tsv, crg2/units.tsv, crg2/config_hpf.yaml, crg2/dnaseq_slurm_hpf.sh, crg2/slurm_profile/slurm-config.yaml .
You may need to re-copy config_hpf.yaml and slurm-config.yaml if the files were recently updated in the repo from previous crg2 runs. Note that 'slurm-config.yaml' is for submitting each rule as cluster jobs, so ignore this if not running on cluster.
```
mkdir NA12878
cp crg2/samples.tsv crg2/units.tsv crg2/config_hpf.yaml crg2/dnaseq_slurm_hpf.sh crg2/slurm_profile/slurm-config.yaml NA12878
cd NA12878
```

2. Set up pipeline run 
  * Reconfigure 'samples.tsv', 'units.tsv' to reflect sample names and input files.
  * Modify 'config_hpf.yaml': 
    * change `project` to refer to the family ID (here, NA12878).
    * set `pipeline` to `wes`, `wgs`, `annot` or `mity` for exome sequences, whole genome sequences, to simply annotate a VCF respectively or to generate mitochondrial reports
    * inclusion of a panel bed file (`panel`) or hpofile (`hpo`) will generate 2 SNV reports with all variants falling within these regions, one which includes variants in flanking regions as specified by `flank`.
    * inclusion of a `ped` file with parents and a proband(s) will allow generation of a genome-wide de novo report if `pipeline` is `wgs`.
    * `minio` refers to the path of the file system mount that backs MinIO in the CHEO HPC4Health tenancy. Exome reports (and coverage reports, if duplication percentage is >20%) will be copied to this path if specified. 

samples.tsv
```
sample
NA12878
```

units.tsv
```
sample	platform	fq1	fq2	bam	cram
NA12878	ILLUMINA			/hpf/largeprojects/ccmbio/GIAB_benchmark_datasets/hg19/WGS/NA12878/RMNISTHS_30xdownsample.bam
```

config_hpf.yaml
```
run:
  project: "NA12878"
  samples: samples.tsv
  units: units.tsv
  ped: "" # leave this string empty if there is no ped
  panel: "" # three-column BED file based on hpo file; leave this string empty if there is no panel
  hpo: "" # five-column TSV with HPO terms; leave this string empty is there are no hpo terms
  flank: 100000
  gatk: "gatk"
  pipeline: "wgs" #either wes (exomes) or wgs (genomes) or annot (to annotate and produce reports for an input vcf) or mity (to generate mitochondrial reports)
  minio: ""
  PT_credentials: ""
...
```

3. Activate the conda environment with Snakemake 5.10.0

```
(base) [dennis.kao@qlogin5 crg2]$ conda activate snakemake
(snakemake) [dennis.kao@qlogin5 crg2]$ snakemake -v
5.10.0
```

4. Test that the pipeline will run by adding the flag "-n" to the command in dnaseq_slurm_hpf.sh and running it.

```
(snakemake) [dennis.kao@qlogin5 crg2]$ sh dnaseq_slurm_hpf.sh
Building DAG of jobs...
Job counts:
	count	jobs
	1	EH
	1	EH_report
	1	EHdn
	1	EHdn_report
	1	all
	1	allsnvreport
	1	bam_to_cram
	1	bcftools_stats
	1	bgzip
	25	call_variants
	25	combine_calls
	1	fastq_screen
	1	fastqc
	25	genotype_variants
	2	hard_filter_calls
	1	input_prep
	1	manta
	1	map_reads
	1	mark_duplicates
	1	md5
	1	merge_calls
	1	merge_variants
	1	metasv
	1	mito_vcfanno
	1	mity_call
	1	mity_normalise
	1	mity_report
	1	mosdepth
	1	multiqc
	1	pass
	1	peddy
	1	qualimap
	1	recalibrate_base_qualities
	1	remove_decoy
	3	samtools_index
	1	samtools_index_cram
	1	samtools_stats
	2	select_calls
	1	smoove
	2	snpeff
	1	subset
	1	svreport
	2	svscore
	1	tabix
	1	vcf2db
	1	vcfanno
	1	vep
	1	verifybamid2
	1	vt
	1	wham
	1	write_version
	129

  [rule inputs and outputs removed for brevity]

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
...
```

5. Submit pipeline job to the cluster

    Job scheduler: PBS on SickKids hpf
      * Serialized jobs:
      ```qsub dnaseq.pbs```
      * Parallelized jobs:
       ```qsub dnaseq_cluster.pbs```

    Refer to `pbs_profile/cluster.md` document for detailed documentation for cluster integration.

    Job scheduler: Slurm
      * Parallelized jobs on SickKids hpf:
      ```sbatch dnaseq_slurm_hpf.sh```
      * Parallelized jobs on CHEO-RI space:
      ```sbatch dnaseq_slurm_cheo_ri.sh```


    `dnaseq_slurm_api.sh` is called by [Stager](https://stager.genomics4rd.ca/) when exome analyses are requested. It automatically sets up (`exome_setup_stager.py`) and kicks off the crg2 WES pipeline via the slurm API using linked files that have been uploaded to MinIO. 

## Pipeline reports

## Reports

Column descriptions and more info on how variants are filtered can be found in the [CCM Sharepoint.](https://sickkidsca.sharepoint.com/:f:/r/sites/thecenterforcomputationalmedicineworkspace/Shared%20Documents/C4Rare-updates/Report_documentation/Documentation?csf=1&web=1&e=nSMt2p)

The WES pipeline generates 5 reports:

1. wes.regular - report on coding SNVs in exonic regions

2. wes.synonymous - report on synonymous SNVs in exonic regions

3. clinical.wes.regular - same as wes.regular but with more stringent filters

4. clinical.wes.synonymous - same as wes.synonymous but with more stringent filters

5. wes.mosaic - putative mosaic variants

The WGS pipeline generates up to 11 reports:

The SNV reports can be found in the directories: 
  - report/coding/{PROJECT_ID}/{PROJECT_ID}.\*wes\*.csv
  - report/panel/{PROJECT_ID}/{PROJECT_ID}.\*wgs\*.csv
  - report/panel-flank/{PROJECT_ID}/{PROJECT_ID}.\*wgs\*.csv
  - report/denovo/{PROJECT_ID}/{PROJECT_ID}.\*wgs\*.csv

The SV reports can be found in the directory (SV, unfiltered SV, and BND): 
  - report/sv/{PROJECT_ID}.wgs.{VER}.{DATE}.tsv

The STR reports can be found in:
  - report/str/{PROJECT_ID}.EH.v1.1.{DATE}.xlsx
  - report/str/{PROJECT_ID}.EHDN.{DATE}.xlsx

The mitochondrial report can be found in the directory:
  - report/mitochondrial/{PROJECT_ID}.mitochondrial.report.csv

The MELT mobile element report can be found in the directory: 
  - report/MELT/{PROJECT_ID}.wgs.{DATE}.csv

## Automatic pipeline submission
### Genomes
`parser_genomes.py` script can be used to automate the above process for a batch of genomes. 

```
usage: parser_genomes.py [-h] -f FILE -s {fastq,mapped,recal,decoy_rm} -d path
Reads sample info from TSV file (-f) and creates directory (-d) necessary to start crg2 from step (-s) requested.
optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Five column TAB-seperated sample info file; template sample file: crg2/sample_info.tsv
  -s {fastq,mapped,recal,decoy_rm}, --step {fastq,mapped,recal,decoy_rm}
                        start running from this folder creation(step)
  -d path, --dir path   Absolute path where crg2 directory struture will be created under familyID as base directory
```
The script performs the following operations for each familyID present in the sample info file
  - create run folder: \<familyID\>
  - copy the following files to run folder and update settings where applicable:
    - config_hpf.yaml: update run name, input type, panel and ped (if avaialble in /hpf/largeprojects/ccmbio/dennis.kao/gene_data/{HPO,Pedigrees})
    - units.tsv: add sample name, and input file paths
    - samples.tsv: add sample names
    - dnaseq_slurm_hpf.sh: rename job 
    - slurm-config.yaml
  - submit Snakemake job 

### Exomes for re-analysis
`exome_reanalysis.py` script can be used to automate the above process for exome re-analysis. 
```
usage: exome_reanalysis.py [-h] -f FILE -d path

Reads sample info from Stager analysis csv file (-f) and creates directory (-d) necessary to run crg2.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Analyses csv output from STAGER
  -d path, --dir path   Absolute path where crg2 directory structure will be created under familyID/analysisID as base directory
  ```

The script parses an analysis request csv from Stager for exome re-analyses and sets up necessary directories (under 2nd argument), files as below:
1. create family analysis directory under directory provided
3. copy config_cheo_ri.yaml, slurm-config.yaml and dnaseq_slurm_cheo_ri.sh from crg2 repo and replace necessary strings
4. search results directory /srv/shared/hpf/exomes/results for crams from previous analyses
5. create units.tsv and samples.tsv for snakemake
6. submit job if all the above goes well

### Genomes for re-analysis
`genome_reanalysis.py` script can be used to automate the above process for genome re-analysis. 
```
usage: genome_reanalysis.py [-h] -f FILE -d path

Reads sample info from Stager analysis csv file (-f) and creates directory (-d) necessary to run crg2.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Analyses csv output from STAGER
  -d path, --dir path   Absolute path where crg2 directory structure will be created under familyID as base directory
  ```

The script parses an analysis request csv from Stager for genome re-analyses and sets up necessary directories (under 2nd argument), files as below:
1. create family analysis directory under directory provided
3. copy config_cheo_ri.yaml, slurm-config.yaml and dnaseq_slurm_cheo_ri.sh from crg2 repo and replace necessary strings
4. search results directories for crams from previous analyses
5. create units.tsv and samples.tsv for snakemake
6. submit job if all the above goes well


## Extra targets

The following output files are not included in the main Snakefile and can be requested in `snakemake` command-line.

1. HPO annotated reports: 
  Reports from coding, panel and panel-flank can be annotated with HPO terms whenever HPO file is available with us. This is done mainly for monthly GenomeRounds. HPO annotated TSV files are created in directory: `report/hpo_annotated` using the following command:

  `snakemake --use-conda -s $SF --conda-prefix $CP --profile ${PBS} -p report/hpo_annotated` 
  
 The reason to not include this output in Snakefile is that the output of `rule allsnvreport` is a directory and hence snakemake will not check for the creation of final csv reports. So, users are required to make sure the csv reports are created in the three folders above, and then request for the output of `rule annotate_hpo` separately on command-line (or append to the dnaseq_cluster.pbs).
  
## Other outputs
CNV and SV comparison outputs are not yet part of the pipeline. Please follow the steps 7 & 8 in [crg](https://github.com/ccmbioinfo/crg#7-cnv-report) to generate the following three TSVs (this is also required for GenomeRounds)
1. \<FAMILYID>.\<DATE>.cnv.withSVoverlaps.tsv
2. \<FAMILYID>.unfiltered.wgs.sv.\<VER>.\<DATE>.withCNVoverlaps.tsv
3. \<FAMILYID>.wgs.sv.\<VER>.\<DATE>.withCNVoverlaps.tsv

## Benchmarking

SNV calls from WES and WGS pipeline can be benchmarked using the GIAB dataset _HG001_NA12878_ (family_sample) and truth calls from NISTv3.3.2

1. Copy all required files for run as [here](#running-the-pipeline).
  The inputs in `units.tsv` are downsampled for testing purposes. Edit the tsv to use the inputs from HPF: `ccmmarvin_shared/validation/benchmarking/benchmark-datasets` or CHEO-RI: `/srv/shared/data/benchmarking-datasets/HG001`. 
2. Copy `crg2/benchmark_hpf.tsv` or `crg2/benchmark_cheo_ri.tsv` to current directory. _Note_: benchmark_x.tsv uses HG001_NA12878 as family_sample name, so you should edit the "project" name in `config_hpf.yaml`
3. Edit the `config_hpf.yaml` to set "wes" or "wgs" pipeline
4. Edit `dnaseq_slurm_hpf.sh` to include the target `validation/HG001`


## crg2 Developer Documentation:

A guideline to developing crg2 can be found in this [document.](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/Shared%20Documents/C4Rare-updates/crg2-development.docx?d=w5cad9def56a44b2bba951e6a0a7a4334&csf=1&web=1&e=FHxe4j)
