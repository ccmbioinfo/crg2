#!/bin/bash

DIR=/hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/914/recal/gene_coverage/

## Row bind gene coverage TSV files for each gene into one file for each sample &
##      add column headers
for sample in "914_CH0528" "914_CH0529" "914_CH0530"; do
    cat ${DIR}"$sample".cram_*.regions.tsv > ${DIR}"$sample".cram_regions.tsv;
    sed  -i '1i CHR\tSTART\tEND\tID\tGENE\tTYPE\tCOVERAGE' ${DIR}"$sample".cram_regions.tsv;
done
