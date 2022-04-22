#!/bin/bash
#SBATCH -A slurm
#SBATCH --partition=all
#SBATCH --job-name=crg2
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out


# sample commands to be called by stager slurm.py with variable substitution
# assumes crg2 is installed in home directory of user
# in production, crg2 is located in /home/ccmmarvin

SF=~/crg2/Snakefile; 
CP="/home/slurm/conda_envs/crg2-conda";
SLURM=~/crg2/slurm_profile;

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate snakemake

family=2910
analysis_id=1235
data_dict=/storage/data/test_crg2_automation/2910.json
filepath=/storage/data/test_crg2_automation
analysis_dir=${filepath}/${analysis_id}/${family}

if [ ! -d $analysis_dir ];then
    mkdir -p $analysis_dir
else
    echo "Analysis directory ${filepath}/${family} exists, exiting"
    exit 1
fi

cd $analysis_dir

python3 ~/crg2/exome_setup_stager.py \
    -a $analysis_dir \
    -f $family \
    -d $data_dict

exit_code=`echo $?`
if [ $exit_code == 0 ]; then
    snakemake --use-conda -s $SF --conda-prefix $CP  --profile $SLURM -p 
else
    echo 'Analysis setup failed, exiting'
    exit
fi