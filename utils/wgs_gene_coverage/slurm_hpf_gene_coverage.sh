#!/bin/bash
#SBATCH --job-name=gene_coverage
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --constraint=CentOS7
#SBATCH --mail-type=ALL

# DIR: path to cram and cram.crai files
# GENE_LIST: list of genes to get coverage for
DIR=/hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/914/recal/
GENE_LIST=/hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/914/recal/gene_coverage/gene_list.txt

## Submit one job for each gene
while read -r line; do 
    echo $line;
    sbatch get_gene_coverage.sh $line $DIR;
done < "$GENE_LIST"
