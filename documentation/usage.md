# Usage

The rule `all` in file `crg2/Snakefile` requests for all the outputs produced by this pipeline, conditional on the settings found in  `config.yaml`. You can requests for specific outputs by passing the exact filenames or directory names as a command-line argument to snakemake. 

## Examples
Do this from the directory you wish to run Snakemake, and check if the dry-run is error free before submitting

0. set up files, variables & activate snakemake env (taken from the `dnaseq_cluster.pbs` file)
   
```bash
mkdir 1155 && cd 1155
cp ~/crg2/config.yaml ~/crg2/samples.tsv ~/crg2/units.tsv ~/crg2/pbs_profile/pbs_config.yaml ~/crg2/danseq_cluster.pbs
sed -i 's/NA12878/1155/' config.yaml #change the project name
SF=~/crg2/Snakefile; 
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake";
PBS=~/crg2/pbs_profile;
source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

```

1. snv coding reports
```bash
snakemake --use-conda -s $SF --conda-prefix $CP  --profile $PBS -p report/coding -n
```

2. for panel reports, set the path of the HPO file in your project 'config.yaml'. (For Care4Rare HPO files needs to be run through a Rscript manually to get a BED file format)
```
run:
  panel: "115R.bed"
```

```bash
snakemake --use-conda -s $SF --conda-prefix $CP  --profile $PBS -p report/panel report/panel-flank -n
```

3. sv reports 
```bash
snakemake --use-conda -s $SF --conda-prefix $CP  --profile $PBS -p report/sv -n
```

4. filtered VCF for validations
```bash
snakemake --use-conda -s $SF --conda-prefix $CP  --profile $PBS -p annotated/coding/vcfanno/all.coding.vep.vcfanno.vcf -n

```

5. multiple files with shell expansion
```bash
snakemake --use-conda -s $SF --conda-prefix $CP  --profile $PBS -p annotated/{coding,panel}/vcfanno/all.{coding,panel}.vep.vcfanno.vcf -n
```
