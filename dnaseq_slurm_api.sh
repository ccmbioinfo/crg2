# sample commands to be called by stager slurm.py with variable substitution
# assumes crg2 is installed in home directory of user
# in production, crg2 is located in /home/ccmmarvin
SF=~/crg2/Snakefile; 
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake";
SLURM=~/crg2/slurm_profile;

source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

if [ ! -d {filepath}/{family} ];then
    mkdir {filepath}/{family}
else
    echo "Analysis directory {filepath}/{family} exists, exiting"
    exit


mkdir {filepath}/{family}
cd {filepath}/{family}

python3 exome_setup.py \
    -a {filepath}/{family} \
    -f {family} \
    -d {participant}:{linked_files} {participant}:{linked_files} {participant}:{linked_files}

exit_code=`echo $?`
if [ $exit_code == 0 ]; then
    snakemake --use-conda -s $SF --conda-prefix $CP  --profile $SLURM -p 
else
    echo 'Analysis setup failed, exiting'
    exit

