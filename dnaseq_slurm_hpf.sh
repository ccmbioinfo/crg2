#!/bin/bash
#SBATCH --job-name=crg2
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out

SF=~/crg2/Snakefile
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake"
SLURM=~/crg2/slurm_profile/
CONFIG=config_hpf.yaml

source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake
#target annotated/coding/vcfanno/NA12878.coding.vep.vcfanno.vcf

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --configfile ${CONFIG} --profile ${SLURM} -p --dry-run 

# annotated/coding/vcfanno/NA12878.coding.vep.vcfanno.vcf
#/hpf/largeprojects/ccmbio/GIAB_benchmark_datasets/hg19/WGS/NA12878
#/hpf/largeprojects/ccmbio/GIAB_benchmark_datasets/hg19/WGS/NA12878/RMNISTHS_30xdownsample.bam