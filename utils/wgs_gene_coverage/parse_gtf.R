library(data.table)
library(rtracklayer)


gtf <- readGFF("/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.gtf")

## Filter out additional chromosomes 
chr_names <- paste("chr", c(seq(1:22), "X", "Y"), sep="")
gtf <- gtf[which(gtf$seqid %in% chr_names),]


## Export transcripts only
transcripts <- gtf[which(gtf$type == "transcript"),]
export(transcripts, "/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.transcripts.gtf")

## Export exons only
exons <- gtf[which(gtf$type == "exon"),]
export(exons, "/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.exons.gtf")
exons$ID <- NA


##### After running get_intron_gft.sh ---
## (uses "bedtools subtract" (transcript-exon) to extract intron coordinates)
introns <- readGFF("/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.introns.gtf")
introns$type <- "intron" 
introns$exon_id <- NA
introns$exon_number <- NA
export(introns, "/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.introns.gtf")


## Merge exons & introns into one BED file
exons_introns_gtf <- rbind(exons, introns)
export(exons_introns_gtf, "/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.exon+intron.gtf")
# exons_introns_gtf <- readGFF("/Users/shaniawu/Desktop/Care4Rare/data/hg19.ncbiRefSeq.exon+intron.gtf")


# Remove "chr" prefix 
exons_introns_gtf$seqid <- gsub("chr", "", exons_introns_gtf$seqid) 


## Extract columns from GTF for BED file ---
# Note: GTF file provides exon transcript and number for exons ("exon_id"), 
#       but only transcript id for introns

bed <- subset(exons_introns_gtf, select = c("seqid", "start", "end", "type", "gene_id", "transcript_id"))

write.table(bed, "/Users/shaniawu/Desktop/Care4Rare/data/hg19_introns_exons.bed", row.names = F)
# also saved to /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/hg19_introns_exons.bed

