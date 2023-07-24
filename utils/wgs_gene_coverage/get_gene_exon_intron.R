library(data.table)
library(rtracklayer)
library(base)

## Read genome intron/exon coordinates (from parse_gft.R)
genome_bed <- fread("/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/hg19_introns_exons.bed")


Get_Gene_Bed <- function(gene_name, dir){
  ## Print input gene_name
  print(paste("Input gene:", gene_name))
  
  ## Check if gene_name input is valid
  valid_genes <- unique(genome_bed$gene_id)
  if (! gene_name %in% valid_genes){
    stop("Invalid gene name")
  }
  
  ## Filter genome_bed for gene regions 
  gene_bed <- genome_bed[which(genome_bed$gene_id == gene_name),]

  ## Concatenate annotation columns into "id" (mosdepth only allows one metadata column)
  gene_bed$id <- paste(gene_bed$transcript_id, gene_bed$gene_id, gene_bed$type, sep=";") 
  gene_bed <- subset(gene_bed, select = -c(transcript_id, gene_id, type))

  ## Make gene_coverage folder in dir and write BED file to folder
  coverage_dir <- paste(dir, "gene_coverage", sep="/")
  dir.create(coverage_dir)

  write.table(gene_bed, paste(coverage_dir, sprintf("%s_pos.bed", gene_name), sep="/"), ######## TODO - check path 
    row.names=F, col.names=F, sep="\t", quote = FALSE)

  print("BED file generated")
}


## Read gene & directory argument from command line
gene <- commandArgs(trailingOnly = TRUE)[2]
project_dir <- commandArgs(trailingOnly = TRUE)[3]

Get_Gene_Bed(gene_name = gene, dir = project_dir)


## Examples
# BRAF_pos <- Get_Gene_Bed('BRAF')
# IRF6_pos <- Get_Gene_Bed('IRF6')
