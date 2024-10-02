rm(list=ls())

#### parameters
#### to run
#### Rscript mergeExpansions.R --rscript ~/crg2/scripts/str/ExpansionAnalysis.R --ehdn /hpf/largeprojects/tcagstor/users/btrost/papers/STRs/Qatar_ASD/output_regions.min2.1000G+SSC+MSSNG+QASD.txt 
#### --outlier dbscan.tsv --outpath output/ --trf /hpf/largeprojects/ccm_dccforge/dccdipg/Common/annotation/ExpansionHunterDenovo/UCSC_simple_repeats_hg19_coord_motif.tsv
###trf <- "/hpf/largeprojects/ccmbio/aarthi/old_projects/proj_CHEO/CRG/str/ExpansionAnalysisPackage/UCSC_simple_repeats_hg19_coord_motif.tsv"
###source(path.expand("~/crg2/str/ExpansionAnalysisFunctions.R"))

#### read arguments from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop(call. = T, "Require more argument(s)")
}

paramNames <- grep("--", args)
paramValues <- paramNames + 1

if(length(grep("--", args[paramValues])) > 0){
  stop(call. = T, "Failed in reading arguments, '--' found in argument value")
}

if(length(args) != length(c(paramNames, paramValues))){
  message(sprintf("Ignore argument(s):%s", paste(args[-union(paramNames, paramValues)], collapse=",")))
}

params <- list()
paramNames <- gsub("--", "", args[paramNames])
for(i in 1:length(paramNames)){
  params[paramNames[i]] <- args[paramValues[i]]
}

source(path.expand(params$rscript))
outpath <- params$outpath
if(length(grep("\\/$", outpath)) == 0){
  outpath <- paste0(outpath, "/")
}

if(!dir.exists(outpath)){
  dir.create(outpath, recursive = T)
}

###############

library(data.table)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(Biostrings)

queryString <- function(str1, str2){
  ori.str2 <- str2
  str2 <- paste0(str2, str2)
  
  length.str1 <- nchar(str1)
  length.str2 <- nchar(str2)
  revcom.str1 <- as.character(reverseComplement(DNAString(str1)))
  
  length.diff <- ifelse(length.str2 > length.str1, length.str2 - length.str1, 1)
  
  identity <- c()
  for(i in 1:length.diff){
    tmp <- substr(str2, i, length.str1+i-1)
    match <- sum(strsplit(tmp, "")[[1]] == strsplit(str1, "")[[1]])
    match.recom <- sum(strsplit(tmp, "")[[1]] == strsplit(revcom.str1, "")[[1]])
    
    match <- min(max(match, match.recom), nchar(ori.str2))
    tmp.identity <- match/length.str1
    
    identity <- c(identity, tmp.identity)
  }
  return(max(identity))
}

#### prepare map TRF
message(sprintf("reading %s##", params$ehdn))
ehdn <- fread(params$ehdn)
trf <- fread(params$trf)
ehdn <- as.data.frame(ehdn)[, c(1:6)]
ehdn$repeatID <- paste(ehdn$V1, ehdn$V2, ehdn$V3, ehdn$V4, sep="#")
names(ehdn) <- c("chr", "start", "end", "motif", "var1", "var2", "repeatID")
if(length(grep("chr", ehdn$chr[1:5])) == 0)
  ehdn$chr <- paste0("chr", ehdn$chr)
trf <- as.data.frame(trf)
trf$repeatID <- paste(trf$V1, trf$V2, trf$V3, trf$V4, sep="#")

ehdn.g <- GRanges(ehdn$chr, IRanges(ehdn$start, ehdn$end), "*")
trf.g <- GRanges(trf$V1, IRanges(trf$V2, trf$V3), "*")

olap <- data.frame(findOverlaps(trf.g, ehdn.g))
olap$trf.motif <- trf$V4[olap$queryHits]
olap$ehdn.motif <- ehdn$motif[olap$subjectHits]
olap$trf.id <- trf$repeatID[olap$queryHits]
olap$ehdn.id <- ehdn$repeatID[olap$subjectHits]

for(i in 1:nrow(olap)){
  trf.identity <- queryString(olap$trf.motif[i], olap$ehdn.motif[i])
  ehdn.identity <- queryString(olap$ehdn.motif[i], olap$trf.motif[i])
  
  olap$reciprocalIdentity[i] <- min(trf.identity, ehdn.identity)
}

olap$length.diff <- nchar(olap$trf.motif) - nchar(olap$ehdn.motif)

olap.f <- olap[order(olap$reciprocalIdentity, decreasing = T), ]
olap.f <- olap.f[!duplicated(olap.f$subjectHits), ]
olap.f <- olap.f[olap.f$reciprocalIdentity > 0.66, ]
olap.f <- aggregate(trf.id ~ ehdn.id + reciprocalIdentity + length.diff, olap.f, paste, collapse = ";")
olap.f <- unique(olap.f)
olap.f <- olap.f[order(olap.f$reciprocalIdentity, decreasing = T), ]

olap.f <- olap.f[!duplicated(olap.f$ehdn.id), ]

write.table(olap.f[, c(4, 1, 2, 3)], sprintf("%smap.TRF.EHdn.0.66.tsv", outpath), sep="\t", row.names=F, quote=F, col.names=T)

#########################################
### merge
outfile <- sprintf("%smerged.ehdn.tsv", outpath)

dt <- fread(params$ehdn)
dt <- data.frame(dt)
ehdn <- dt
dt$varid <- paste(dt$V1, dt$V2, dt$V3, sep="_")

clean.sample <- readLines(sprintf("%sclean.samples.txt", outpath))

ehdn$repeatID <- paste(ehdn$V1, ehdn$V2, ehdn$V3, ehdn$V4, sep="#")
ehdn$samples <- 0
for(i in 1:nrow(ehdn)){
  tmp <- strsplit(gsub(",", ":", ehdn$V7[i]), ":")[[1]]
  tmp <- tmp[seq(1, length(tmp), 2)]
  tmp <- tmp[tmp %in% clean.sample]
  
  ehdn$samples[i] <- length(tmp)
}

ehdn <- ehdn[ehdn$samples > 0, ]
map <- read.delim(sprintf("%smap.TRF.EHdn.0.66.tsv", outpath), stringsAsFactors = F)

ehdn.all <- merge(ehdn[, -c(5:7, 9)], map, by.x = "repeatID", by.y = "ehdn.id", all.x = T)
names(ehdn.all) <- c("repeatID", "chr", "start", "end", "motif", "trf.id", "reciprocalIdentity", "length.diff")
ehdn.all <- ehdn.all[ehdn.all$chr %in% paste0("chr", c(1:22, "X", "Y")), ]

ehdn.m <- mergeLoci(ehdn.all)
ehdn.m$uniqueMotif <- sapply(sapply(ehdn.m$motif, strsplit, ";"), length)

ehdn.m$allRepeat <- ehdn.m$repeatID
ehdn.m$repeatID <- paste(ehdn.m$chr, ehdn.m$start, ehdn.m$end, sep="#")

write.table(ehdn.m, outfile, sep="\t", row.names=F, quote=F, col.names=T)

##################
### extract rare merged expansion
expansion <- read.delim(params$outlier, stringsAsFactors = F)
expansion <- expansion[expansion$a1000g_freq_perc < 0.1, ]

ehdn.m.g <- GRanges(ehdn.m$chr, IRanges(ehdn.m$start, ehdn.m$end), "*")
expansion.g <- GRanges(expansion$chr, IRanges(expansion$start, expansion$end), "*")

olap <- data.frame(findOverlaps(expansion.g, ehdn.m.g))
olap <- aggregate(queryHits ~ subjectHits, olap, paste, collapse = ";")

expansion.m <- data.frame()
for(i in 1:nrow(olap)){
  idx <- as.numeric(strsplit(olap$queryHits[i], ";")[[1]])
  tmp.exp <- expansion[idx, ]
  motif <- paste(tmp.exp$motif, sep=";")
  outliers <- paste(tmp.exp$outliers, sep=";")
  repeats <- paste(tmp.exp$repeatID, sep=";")
  chr <- tmp.exp$chr[1]
  start <- min(tmp.exp$start)
  end <- max(tmp.exp$end)
  
  expansion.m <- rbind(expansion.m, data.frame(chr, start, end, motif, outliers, stringsAsFactors = F))
}

write.table(expansion.m, sprintf("%smerged.rare.expansions.tsv", outpath), sep="\t", row.names=F, col.names=T, quote=F)
            
expansion.m$ref <- "0"
expansion.m$alt <- "-"
expansion.m$varid <- paste(expansion.m$chr, expansion.m$start, expansion.m$end, sep="#")
names(expansion.m)[1] <- "#chr"
exp.annovar <- expansion.m[, c("#chr", "start", "end", "ref", "alt", "varid")]

write.table(exp.annovar, sprintf("%smerged.rare.expansions.forannotation.tsv",outpath), sep="\t", row.names=F, quote=F, col.names=T)
#####
#####
###print("Printing sessioninfo from mergeExpansions.R")
###sessionInfo()
