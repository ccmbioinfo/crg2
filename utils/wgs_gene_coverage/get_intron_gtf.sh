DATA_DIR=~/crg2/utils/wgs_gene_coverage/data;

bedtools subtract -a ${DATA_DIR}/hg19.ncbiRefSeq.transcripts.gtf -b ${DATA_DIR}/hg19.ncbiRefSeq.exons.gtf > ${DATA_DIR}/hg19.ncbiRefSeq.introns.gtf