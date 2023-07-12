library(data.table)
library(rtracklayer)

dir <- "~/crg2/utils/wgs_gene_coverage/"


print("Reading GTF file...")
gtf <- readGFF("/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/hg19.ncbiRefSeq.gtf")

## Filter out additional chromosomes 
chr_names <- paste("chr", c(seq(1:22), "X", "Y"), sep="")
gtf <- gtf[which(gtf$seqid %in% chr_names),]

## Make "data" folder
dir.create(paste(dir, "data", sep=""))


## Export transcripts only ---
print("Getting transcripts...")
transcripts <- gtf[which(gtf$type == "transcript"),]
export(transcripts, paste(dir, "data/hg19.ncbiRefSeq.transcripts.gtf", sep="")) # creates

## Export exons only ---
print("Getting exons...")
exons <- gtf[which(gtf$type == "exon"),]
export(exons, paste(dir, "data/hg19.ncbiRefSeq.exons.gtf", sep=""))
exons$ID <- NA


## Get intron coordinates using bedtools subtract (transcript-exon) ---
print("Getting introns...")
system("DATA_DIR=~/crg2/utils/wgs_gene_coverage/data; bedtools subtract -a ${DATA_DIR}/hg19.ncbiRefSeq.transcripts.gtf -b ${DATA_DIR}/hg19.ncbiRefSeq.exons.gtf > ${DATA_DIR}/hg19.ncbiRefSeq.introns.gtf")

introns <- readGFF(paste(dir, "data/hg19.ncbiRefSeq.introns.gtf", sep=""))
introns$type <- "intron" 
introns$exon_id <- NA
introns$exon_number <- NA


## Merge exons & introns into one BED file & format ---
print("Formatting and creating exons+introns BED...")
exons_introns_gtf <- rbind(exons, introns)
exons_introns_gtf$seqid <- gsub("chr", "", exons_introns_gtf$seqid) # Remove "chr" prefix 

bed <- subset(exons_introns_gtf, select = c("seqid", "start", "end", "type", "gene_id", "transcript_id"))


## Write BED file ---
write.table(bed, (paste(dir, "data/hg19_introns_exons.bed", sep="")), row.names = F)
# copy of file in /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/hg19_introns_exons.bed used in get_gene_exon_intron.R

print("Done. BED file saved to hg19_introns_exons.bed.")


## Delete interediate files
del <- c(paste(dir, "data/hg19.ncbiRefSeq.transcripts.gtf", sep=""),
        paste(dir, "data/hg19.ncbiRefSeq.exons.gtf", sep=""),
        paste(dir, "data/hg19.ncbiRefSeq.introns.gtf", sep=""))
unlink(del)
print("Deleted intermediate files.")