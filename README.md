# crg2-bacterial
Research pipeline for de-novo assembly of whole genome (WGS) of bacterial isolates. It aims to provide reproducible results, be computationally efficient, and transparent in it's workflow.crg2 uses Snakemake and Conda to manage jobs and software dependencies.

<div align="center">
    <img src="/crg2logolarge.png" width="800px"</img> 
</div>

## Conda environment set-up
Snakemake will build the conda-environment for each rule according to the environment.yaml file present in the respective wrapper's directory. For some tools, we need to download the database separately.

### Setting up the Kraken2 database
1. Download pre-built database: `wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250714.tar.gz`
2. Create a directory: `mkdir /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/k2_standard`
3. Unzip the file: `tar -xvf k2_standard_20250714.tar.gz -C k2_standard`
4. Add the database path to the `config_hpf.yaml`

### Setting up the Krona database 
1. Activate conda environment created by snakemake: `conda activate /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/0017e84a`
2. Remove the taxonomy directory:  `rm -rf /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/0017e84a/opt/krona/taxonomy`
3. Create the taxonomy directory in another folder:  `mkdir /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/krona/taxonomy`
4. Create a symlink: `ln -s /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/krona/taxonomy /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/0017e84a/opt/krona/taxonomy`
5. Change directory to the newly created taxonomy directory: `cd /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/krona/taxonomy`
6. Download the taxonomy file: `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
7. Build the taxonomy database:  `ktUpdateTaxonomy.sh --only-build`

### Setting up the Bakta database
To set up the Bakta environment, I had to follow a few additional steps. Documenting it here for future reference:
1. Activate the conda environment that snakemake created: `conda activate /hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake/26ae42d8`
2. Change directory to: `cd /hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/data/`
3. Download the bakta database:  `wget https://zenodo.org/records/14916843/files/db.tar.xz?download=1`
4. Rename the downloaded file: `mv db.tar.xz?download=1 db.tar.xz`
5. Run: `bakta_db install -i db.tar.xz`
6. Add the database path to the `config.yaml`


## Pipeline details

### Gather input file
1. Download the run file associated to the accession provided in `accession.tsv` using the SRA toolkit (prefetch)

2. Extract fasta files from SRA run files using the SRA toolkit (fasterq-dump).

### Pre-assembly QC
1. Perform QC using fastqc, kraken2, krona and bracken2 on the raw fastq reads.

### Assembly
1. Perform assembly using Shovill (spades assembler)

### Post-assembly QC
1. Perform QC using kraken2, krona and bracken2 on the shovill assembled assembly.

### Annotation
1. Annotate the assembly using Bakta


## Running the pipeline
1. Make a folder in a directory with sufficient space. Copy over the template files crg2/accession.tsv, crg2/config_hpf.yaml, crg2/dnaseq_slurm_hpf.sh, crg2/slurm_profile/slurm-config.yaml .
You may need to re-copy config_hpf.yaml and slurm-config.yaml if the files were recently updated in the repo from previous crg2 runs. Note that 'slurm-config.yaml' is for submitting each rule as cluster jobs, so ignore this if not running on cluster.
```
mkdir SRS5146886
cp crg2/accession.tsv crg2/config_hpf.yaml crg2/dnaseq_slurm_hpf.sh crg2/slurm_profile/slurm-config.yaml SRS5146886
cd SRS5146886
```

2. Set up pipeline run 
  * Reconfigure 'accession.tsv' to reflect SRA run accession.
  * Modify 'config_hpf.yaml': 
    * change `project` to refer to the SRA sample (here, SRS5146886).

accession.tsv
```
sra_run
SRR9824559
```

config_hpf.yaml
```
run:
  project: "SRS5146886"
  accession: accession.tsv
...
```

3. Activate the conda environment with Snakemake 5.10.0

```
(base) [ajain@login4 SRS5146886]$ source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
(base) [ajain@login4 SRS5146886]$ conda activate snakemake
(snakemake) [ajain@login4 SRS5146886]$ snakemake -v
5.10.0

```

4. Test that the pipeline will run by adding the flag "-n" to the command in dnaseq_slurm_hpf.sh and running it.

```
(snakemake) [dennis.kao@qlogin5 crg2]$ sh dnaseq_slurm_hpf.sh
Building DAG of jobs...
Checking status of 0 jobs.
Using shell: /usr/bin/bash
Provided cluster nodes: 4
Provided resources: cpus=64
Job counts:
        count   jobs
        1       all
        1       bakta
        1       bracken_post_assembly
        1       bracken_pre_assembly
        1       fasterq_dump
        1       fastqc
        1       kraken2_post_assembly
        1       kraken2_raw_reads
        1       krona_post_assembly
        1       krona_pre_assembly
        1       prefetch_sra
        1       quast
        1       shovill
        13


  [rule inputs and outputs removed for brevity]

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
...
```

5. Submit pipeline job to the cluster

    Job scheduler: Slurm
      * Parallelized jobs on SickKids hpf:
      ```sbatch dnaseq_slurm_hpf.sh```

# References:
1. SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump

2. Kraken2: https://github.com/DerrickWood/kraken2

3. Krona: https://github.com/marbl/Krona/wiki

4. Bracken: https://github.com/jenniferlu717/Bracken

5. Shovill: https://github.com/tseemann/shovill

6. Bakta: https://github.com/oschwengers/bakta
