#!/bin/bash
#SBATCH --job-name=gene_coverage
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --constraint=CentOS7
#SBATCH --mail-type=ALL

set -e 

GENE=$1
DIR=$2 # DIR to folder containing cram files to be analyzed 

# Run R script to retrieve exon and intron positions 
module load R
Rscript ./get_gene_exon_intron.R --args ${GENE} ${DIR}


# Run mosdepth ---
export REF_PATH=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa

# $i = cram file full name 
# $prefix = cram file name (only family_sample)
for i in `ls ${DIR}/*.cram`; do
    # Run mosdepth for each cram file in folder
    prefix=`echo $i | cut -d'-' -f1`; 
    mosdepth --by ${DIR}/gene_coverage/${GENE}_pos.bed \
            --thresholds 10,20 \
            "${prefix}_${GENE}" \
            $i;
    echo "Mosdepth running for ${prefix}_${GENE}..."
done;


## Unzip *.cram.regions.bed.gz file 
gunzip ${DIR}*_${GENE}.regions.bed.gz

# Move file into ./gene_coverage/ folder 
mv "${DIR}"*_"${GENE}".regions.bed ${DIR}/gene_coverage/;

# Delete additional files in DIR
rm -r "${DIR}"*_${GENE}.mosdepth.*.dist.txt "${DIR}"*_${GENE}.per-base.bed.gz* "${DIR}"*_${GENE}.regions.bed.gz.csi "${DIR}"*_${GENE}.thresholds.bed.gz* 


## Format by removing ";" and changing to column (SAMPLE_FAM.GENE.cram.regions.bed)
OUT_BED="${DIR}"gene_coverage/*_"${GENE}".regions.bed

for file in ${OUT_BED}; do
    # format columns
    awk '{gsub (/;/, OFS)} 1' OFS="\t" $file > $file.tmp && mv $file.tmp $file;
    # change extension
    tsv=$(echo $file | sed -e 's/bed/tsv/g');
    mv $file $tsv;
done


## Check if mosdepth output created & format (SAMPLE_FAM.GENE.cram.regions.bed)
OUT_TSV=${DIR}/gene_coverage/*_${GENE}.regions.tsv

for file in ${OUT_TSV}; do
    if [[ -f $file ]]; then
        echo "Mosdepth coverage complete, saved to $file";
    fi
done
