#!/usr/bin/env bash
# assumes crg2 is installed in home directory of user
# in production, crg2 is located in /home/ccmmarvin
# ANALYSIS_ID, FAMILY, and DATA are provided by STAGER

set -euo pipefail

SF=~/crg2/Snakefile; 
CP="/home/slurm/conda_envs/crg2-conda";
SLURM=~/crg2/slurm_profile;
ANALYSIS_ID="$1"
FAMILY="$2"
DATA="$3"
#temporary filepath
FILEPATH=/storage/data/test_crg2_automation
ANALYSIS_DIR=${FILEPATH}/${ANALYSIS_ID}/${FAMILY}

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate snakemake

if [ ! -d $ANALYSIS_DIR ];then
    mkdir  -p $ANALYSIS_DIR
else
    echo "Analysis directory $ANALYSIS_DIR exists, exiting"
    exit 1
fi

cd $ANALYSIS_DIR

python3 exome_setup.py \
    -a $ANALYSIS_DIR \
    -f $FAMILY \
    -d $DATA

exit_code=`echo $?`
if [ $exit_code == 0 ]; then
    snakemake --use-conda -s $SF --conda-prefix $CP  --profile $SLURM -p 
else
    echo 'Analysis setup failed, exiting'
    exit 1
fi
