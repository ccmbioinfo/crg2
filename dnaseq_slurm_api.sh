# sample commands to be called by stager slurm.py
# maybe 
SF=~/crg2/Snakefile; 
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake";
SLURM=~/crg2/slurm_profile;

source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

mkdir {filepath}/{family}
cd {filepath}/{family}

python3 exome_setup.py \
    -a {filepath}/{family} \
    -f {family} \
    -d {participant}:{linked_files} {participant}:{linked_files} {participant}:{linked_files}

snakemake --use-conda -s $SF --conda-prefix $CP  --profile $SLURM -p 