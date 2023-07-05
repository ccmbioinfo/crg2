#!/bin/bash
#SBATCH --job-name=gene_coverage
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --constraint=CentOS7
#SBATCH --mail-type=ALL

set -e 

GENE=$1
DIR=$2 # DIR to folder containing cram files to be analyzed 
SCRIPT_DIR=~/crg2/utils/wgs_gene_coverage 
# Run R script to retrieve exon and intron positions 
module load R
Rscript ${SCRIPT_DIR}/get_gene_exon_intron.R --args ${GENE} ${DIR}


# Run mosdepth ---
export REF_PATH=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa

# $i = cram file full name 
# $prefix = cram file name (only family_sample)
for i in `ls ${DIR}/*.bam`; do
    echo $i
    # Run mosdepth for each cram file in folder
    prefix=`echo "${i/.bam/}"`
    echo "Mosdepth running for ${prefix}_${GENE}..."
    mosdepth --by ${DIR}/gene_coverage/${GENE}_pos.bed \
            "${prefix}_${GENE}" \
            $i;
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
