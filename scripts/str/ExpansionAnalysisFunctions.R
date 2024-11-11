get1000Freq <- function(outlier, fam.1000g){
  outlier <- strsplit(outlier, ";")[[1]]
  count <- sum(unique(outlier) %in% fam.1000g)
  return(round(count/length(fam.1000g) * 100, digits = 2))
}


mergeLoci <- function(loci, mergeAny = T){
  
  cnv.tmp <- loci
  
  firstIt <- T
  olap <- data.frame()
  while(nrow(olap) != 0 | firstIt){
    firstIt <- F
    cnv.query <- GRanges(cnv.tmp$chr, IRanges(cnv.tmp$start, cnv.tmp$end), "*")
    
    olap <- data.frame(findOverlaps(cnv.query, cnv.query))
    olap <- olap[which(olap$queryHits < olap$subjectHits), ]
    olap <- olap[!duplicated(olap$queryHits), ]
    olap <- olap[!olap$queryHits %in% olap$subjectHits, ]
    if(!mergeAny){
      olap <- olap[nchar(cnv.tmp$motif[olap$queryHits]) == nchar(cnv.tmp$motif[olap$subjectHits]), ]
    }
    
    message(sprintf("%s found overlapped, from %s loci", nrow(olap), nrow(cnv.tmp)))
    flush.console()
    
    if(nrow(olap) > 0){
      ids <- union(olap$queryHits, olap$subjectHits)
      merge.dt <- data.frame()
      for(i in 1:nrow(olap)){
        olap.th <- olap[i, ]
        
        olap.th[, c("qstart", "qend")] <- cnv.tmp[olap.th$queryHits, c("start", "end")]
        olap.th[, c("sstart", "send")] <- cnv.tmp[olap.th$subjectHits, c("start", "end")]
        
        
        repeatID <- paste(union(strsplit(cnv.tmp$repeatID[olap.th$queryHits], ";")[[1]], 
                                strsplit(cnv.tmp$repeatID[olap.th$subjectHits], ";")[[1]]), collapse=";")
        motif <- paste(union(strsplit(cnv.tmp$motif[olap.th$queryHits], ";")[[1]], 
                             strsplit(cnv.tmp$motif[olap.th$subjectHits], ";")[[1]]), collapse = ";")
        trf.id <- paste(na.omit(union(strsplit(cnv.tmp$trf.id[olap.th$queryHits], ";")[[1]], 
                                      strsplit(cnv.tmp$trf.id[olap.th$subjectHits], ";")[[1]])), collapse = ";")
        reciprocalIdentity <- max(cnv.tmp$reciprocalIdentity[olap.th$queryHits], cnv.tmp$reciprocalIdentity[olap.th$subjectHits], na.rm =T)
        length.diff <- 0
        
        merge.dt <- rbind(merge.dt, data.frame(repeatID, "chr" = cnv.tmp$chr[olap.th$queryHits],
                                               "start" = pmin(olap.th$qstart, olap.th$sstart), "end" = pmax(olap.th$qend, olap.th$send),
                                               motif, trf.id, reciprocalIdentity, length.diff))
      }
      
      
      cnv.tmp <- cnv.tmp[-(ids), ]
      cnv.tmp <- rbind(cnv.tmp, merge.dt)
    }
    
    gc();
  }
  
  return(cnv.tmp)
}

mapOutlierToMergeLoci <- function(outlier.loci, merge.loci){
  tmp <- unique(as.numeric(unlist(sapply(outlier.loci$repeatID, grep, merge.loci$repeatID))))
  merge.loci <- merge.loci[tmp, ]
  
  for(i in 1:nrow(merge.loci)){
    repeatIDs <- strsplit(merge.loci$repeatID[i], ";")[[1]]
    outliers <- paste(unique(unlist(sapply(outlier.loci$outliers[outlier.loci$repeatID %in% repeatIDs], strsplit, ";"))), collapse = ";")
    merge.loci$outliers[i] <- outliers
  }
  
  return(merge.loci)
}

getOutlierMergeData <- function(i, ehdn.result.m, fam.data){
  ehdn.rec <- ehdn.result.m[i, ]
  outliers <- strsplit(ehdn.rec$outliers, ";")[[1]]
  
  
  ehdn.data <-  fam.data[fam.data$ID %in% outliers, c("ID", "Affection")]
  if(nrow(ehdn.data) == 0){
    return(NULL)
  }else{
    
    ehdn.data[, c("repeatID", "motif", "chr", "start", "end", "trf.id", "gene_symbol", "entrez_id", "typeseq_priority")] <- 
      ehdn.rec[, c("repeatID", "motif", "chr", "start", "end", "trf.id", "gene_symbol", "entrez_id", "typeseq_priority")]
  }
  return(ehdn.data)
}


getCountTable <- function(dt, label, label.list, fam.data, ehdn.genotype.data, cov){
  if(nrow(dt) == 0){
    return(NA)
  }else{
    dt <- dt[dt$ID %in% fam.data$ID, ]
    
    if(sum(names(fam.data) %in% cov) > 0)
      cov.dt <- unique(fam.data[, c("ID", names(fam.data)[names(fam.data) %in% cov])])
    
    table.tmp <- aggregate(repeatID ~ ID, dt, length)
    names(table.tmp)[2] <- "ExpansionCount" 
    
    table.tmp <- merge(table.tmp, fam.data[, c("ID", label)], by.x = "ID", by.y = "ID", all = T)
    table.tmp[is.na(table.tmp)] <- 0
    table.tmp[, label] <- ifelse(table.tmp[, label] == label.list[1], 1, 0)
    table.tmp[, label] <- factor(table.tmp[, label])
    if(sum(names(fam.data) %in% cov) > 0)
      table.tmp <- merge(table.tmp, cov.dt, by.x = "ID", by.y = "ID", all.x = T)
    table.tmp <- merge(table.tmp, ehdn.genotype.data, by = "ID", all.x = T)
    table.tmp[is.na(table.tmp)] <- 0
    return(table.tmp)
  }
}

testGLM <- function(dt, label, cov, feature, standardization = T){
  if(is.null(nrow(dt))){
    return(list("OR" = NA, "upper" = NA, "lower" = NA,  "pvalue" = NA))
  }else if (length(unique(na.omit(dt[, label]))) < 2){
    return(list("OR" = NA, "upper" = NA, "lower" = NA,  "pvalue" = NA))
  }else{
    for(cov.each in cov){
      if(length(unique(dt[, cov.each])) < 2){
        cov <- cov[-which(cov == cov.each)]
      }
    }
    
    if(standardization){
      dt[, feature] <- scale(dt[, feature])
    }
    ref <- paste0(label, " ~ ", paste(cov, collapse = " + "))
    add <- paste0(ref, " + ", feature)
    ref.lm <- glm(ref, dt, family = binomial(link = "logit"))
    add.lm <- glm(add, dt, family = binomial(link = "logit"))
    ano <- anova(ref.lm, add.lm, test = "Chisq")
    
    or <- exp(add.lm$coefficients[feature])
    pvalue <- ano$`Pr(>Chi)`[2]
    
    conf <- confint(add.lm)
    low <- exp(conf[feature, 1])
    up <- exp(conf[feature, 2])
    
    return(list("OR" = or, "upper" = up, "lower" = low, "pvalue" = pvalue))
  }
}



testByAggregate <- function(repeatIDs, tmp.f, fam.data, ehdn.genotype.data, cov, standardization = F){
  caseCtrl.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & 
                               tmp.f$Affection %in% 1:2), ]
  
  feature <- "ExpansionCount"
  
  all.table <- getCountTable(caseCtrl.dt, "Affection", 
                             c(2, 1), 
                             fam.data, ehdn.genotype.data, cov)
  
  res <- testGLM(all.table, "Affection", cov, feature, standardization)
  
  return(data.frame(res, "affected.count" = sum(all.table$Affection == 1 & all.table$ExpansionCount > 0),
                    "unaffected.count" = sum(all.table$Affection == 0 & all.table$ExpansionCount > 0)))
}
