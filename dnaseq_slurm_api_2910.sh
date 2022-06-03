#!/bin/bash
#SBATCH -A slurm
#SBATCH --partition=all
#SBATCH --job-name=crg2
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out


# sample commands to be called by stager slurm.py with variable substitution

#set -euo pipefail
# comment out above because leads to /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh: line 55: PS1: unbound variable 

SF="/srv/shared/pipelines/crg2/Snakefile"; 
CP="/srv/shared/conda_envs/crg2-conda/";
SLURM="/srv/shared/pipelines/crg2/slurm_profile";

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate /srv/shared/conda_envs/snakemake_5.10.0/

family=2910
analysis_id=1236
data_dict='{"SK0418": ["skh/2190_SK0418/2020-163-144-NGSCUSTOM_markdup.bam","skh/2190_SK0418/2020-163-144-NGSCUSTOM_markdup.bam.bai","skh/2190_SK0418/2020-163-144-NGSCUSTOM_S16_L001_R1_001_211F.fastq.gz","skh/2190_SK0418/2020-163-144-NGSCUSTOM_S16_L001_R2_001_211F.fastq.gz","skh/2190_SK0418/2020-163-144-NGSCUSTOM_S16_L002_R1_001_211F.fastq.gz","skh/2190_SK0418/2020-163-144-NGSCUSTOM_S16_L002_R2_001_211F.fastq.gz"]}'
filepath=/srv/shared/data/test_crg2_automation
analysis_dir=${filepath}/${family}/${analysis_id}

if [ ! -d $analysis_dir ];then
    mkdir -p $analysis_dir
else
    echo "Analysis directory ${analysis_dir} exists"
    exit 1
fi

cd $analysis_dir

python3 ~/crg2/exome_setup_stager.py \
    -a "$analysis_dir" \
    -f "$family" \
    -d "$data_dict"

snakemake --use-conda -s $SF --conda-prefix $CP  --profile $SLURM -p 
