#!/bin/bash

#PBS -N ExpansionHunter
#PBS -l vmem=20g,mem=20g,nodes=1:ppn=1,walltime=24:00:00
#PBS -joe
#PBS -d .

#usage: sh run_str.sh <familyid>
#usage: qsub run_str.sh -F "<familyid>"

family=$1

module load ExpansionHunter/3.0.1

ref=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37d5/seq/GRCh37d5.fa
catalog=/hpf/largeprojects/ccmbio/arun/C4Rare/Tandem_repeat_disease_loci.hg19.json

for reads in ${family}/bcbio-align/${family}/final/${family}_*/${family}_*ready.bam; do 

    if [ ! -d "${family}/str/expansion_hunter" ]; then 
        mkdir -p "${family}/str/expansion_hunter";
    fi;

    if [ ! -d "${family}/str/expansion_hunter_denovo" ]; then 
        mkdir -p "${family}/str/expansion_hunter_denovo";
    fi;
    
    
    prefix=`basename $reads -ready.bam`;
    eh_prefix="${family}/str/expansion_hunter/${prefix}";
    ehdn_prefix="${family}/str/expansion_hunter_denovo/${prefix}";

    #eh
    sex=`sh ~/crg2/scripts/str_helper.sh $reads`;
    echo  -e "ExpansionHunter --reads $reads --reference $ref --variant-catalog $catalog --output-prefix $eh_prefix --sex $sex"
    ExpansionHunter --reads $reads --reference $ref --variant-catalog $catalog --output-prefix $eh_prefix --sex $sex

    #ehdn
    echo -e "/hpf/largeprojects/ccmbio/arun/Tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0 --reference $ref --reads $reads --output-prefix $ehdn_prefix --min-anchor-mapq 50 --max-irr-mapq 40"
    /hpf/largeprojects/ccmbio/arun/Tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0 --reference $ref --reads $reads --output-prefix $ehdn_prefix --min-anchor-mapq 50 --max-irr-mapq 40

done; 
