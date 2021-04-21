#!/bin/bash
#SBATCH -A slurm
#SBATCH --partition=all
#SBATCH --job-name=crg2
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out

SF="/home/slurm/crg2/Snakefile"
CP="/home/slurm/conda_envs/crg2-conda/"
SLURM="/home/slurm/crg2/slurm_profile/"

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate snakemake_5.10.0 

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --profile $SLURM  
