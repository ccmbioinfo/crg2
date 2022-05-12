#!/usr/bin/env bash
# assumes crg2 is installed in home directory of user
# in production, crg2 is located in /home/ccmmarvin
# ANALYSIS_ID, FAMILY, and DATA are provided by STAGER

set -euo pipefail

# set environment variables
export ANNOTSV=/srv/shared/data/AnnotSV_2.1
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# set pipeline variables
SF=/srv/shared/pipelines/crg2/Snakefile; 
CP="/home/slurm/conda_envs/crg2-conda";
SLURM=/srv/shared/pipelines/crg2/slurm_profile;
ANALYSIS_ID="$1"
FAMILY="$2"
# DATA is a json string in the format { participant1: [linked files], participant2: ['linked files']...}
DATA="$3"
#temporary filepath
FILEPATH=/srv/shared/analyses/exomes
ANALYSIS_DIR="${FILEPATH}/${ANALYSIS_ID}/${FAMILY}"

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate /srv/shared/conda_envs/snakemake_5.10.0/

if [ ! -d "$ANALYSIS_DIR" ];then
    mkdir  -p "$ANALYSIS_DIR"
    python3 /srv/shared/pipelines/crg2/exome_setup_stager.py \
        -a "$ANALYSIS_DIR" \
        -f "$FAMILY" \
        -d "$DATA"
fi

cd "$ANALYSIS_DIR"

snakemake --use-conda -s "$SF" --conda-prefix "$CP"  --profile "$SLURM" -p 