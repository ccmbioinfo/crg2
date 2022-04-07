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

family=2621
analysis_id=1234
data_dict="'CH2500':['cheo/Exomes/2621_CH2500/2297645.cram','cheo/Exomes/2621_CH2500/2297645.cram.crai','cheo/Exomes/2297645/2621_CH2500.cram.md5','cheo/Exomes/2621_CH2500/2297645.vcf.gz'] 'CH2501':['cheo/Exomes/2621_CH2501/2297669.cram','cheo/Exomes/2621_CH2501/2297669.cram.crai','cheo/Exomes/2621_CH2501/2297669.cram.md5','cheo/Exomes/2621_CH2501/2297669.vcf.gz'] 'CH2502':['cheo/Exomes/2621_CH2502/2297682.cram','cheo/Exomes/2621_CH2502/2297682.cram.crai','cheo/Exomes/2621_CH2502/2297682.cram.md5','cheo/Exomes/2621_CH2502/2297682.vcf.gz']"
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