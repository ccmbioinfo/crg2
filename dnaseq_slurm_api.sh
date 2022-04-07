# sample commands to be called by stager slurm.py with variable substitution
# assumes crg2 is installed in home directory of user
# in production, crg2 is located in /home/ccmmarvin
# {family}, {analysis_id}, {participant}, {linked_files} are passed to the slurm API by slurm.py

SF=~/crg2/Snakefile; 
CP="/home/slurm/conda_envs/crg2-conda";
SLURM=~/crg2/slurm_profile;

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate snakemake

#temporary filepath
filepath=/storage/data/test_crg2_automation
analysis_dir=${filepath}/{analysis_id}/{family}

if [ ! -d $analysis_dir ];then
    mkdir  -p $analysis_dir
else
    echo "Analysis directory $analysis_dir exists, exiting"
    exit 1
fi

cd $analysis_dir

python3 exome_setup.py \
    -a analysis_dir \
    -f {family} \
    -d {participant}:{linked_files} {participant}:{linked_files} {participant}:{linked_files}

exit_code=`echo $?`
if [ $exit_code == 0 ]; then
    snakemake --use-conda -s $SF --conda-prefix $CP  --profile $SLURM -p 
else
    echo 'Analysis setup failed, exiting'
    exit
fi