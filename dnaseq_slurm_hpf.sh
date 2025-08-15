#!/bin/bash
#SBATCH --job-name=crg2
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out

SF=/hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/crg2/Snakefile
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake"
SLURM=/hpf/largeprojects/ccmbio/ajain/isaac_chantel_project/pipelines/crg2/slurm_profile/
CONFIG=config_hpf.yaml

module purge

source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --configfile ${CONFIG} --profile ${SLURM}
