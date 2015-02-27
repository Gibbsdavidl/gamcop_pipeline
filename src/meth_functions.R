
# Meth-specific functions #


makeUniverse <- function()
{
  require("FDb.InfiniumMethylation.hg19")
  `1471.2164.15.51.s2` <- read.csv("/Volumes/StorageDisk/Meth_DF5/pipeline/basis/1471-2164-15-51-s2.csv", stringsAsFactors=FALSE)
  probenames <- `1471.2164.15.51.s2`$probe[`1471.2164.15.51.s2`[,3] == "keep"]
  hm450 <- get450k()
  probes <- hm450[probenames]
  genes <- getNearestGene(probes)  
  universe <- unique(genes$nearestGeneSymbol)
  universe
}


removeDups <- function(batchMat)
{
  # The batches are in date order
  # so keeping the last one of duplicated
  # samples is done.
  
  ms <- batchMat$Specimen
  ms <- str_sub(ms, 1,9)  
  idx <- which(duplicated(ms))
  ids <- ms[idx]
  toRem <- c()
  for (i in 1:length(idx)) {
    jdx <- which(ms == ids[i])
    if (!str_detect(pattern="M-101-", ids[i])) {
      toRem <- c(toRem,jdx)
    } else if (length(jdx) > 1) {
      toRem <- c(toRem,jdx[-jdx[-length(jdx)]])
    }
  }
  batchMat[-toRem,]
}


methAddBatch <- function(clinMat, batchMat, label)
{
  # take the batch table, and add a row to the clinMat with 
  # rowname label
    
  require(stringr)
  batchMat <- removeDups(batchMat)
  ns <- str_sub(colnames(clinMat), 1,7)
  ms <- batchMat$Specimen
  ms <- str_sub(ms, 3,9)  

  clinMat <- rbind(clinMat, rep(NA, ncol(clinMat)))
  rownames(clinMat)[nrow(clinMat)] <- label
  
  j <- nrow(clinMat)
  for (i in 1:length(ns)) {
    ni <- ns[i]
    if (ni %in% ms) {
      mix <- which(ms == ni)
      clinMat[j,i] <- batchMat$BATCH[mix]
    }
  }
}


mapToGenes <- function(deTable, pcut=0.05)
{
  require(stringr)
  require("FDb.InfiniumMethylation.hg19")
  require(org.Hs.eg.db)
  hm450 <- get450k()
  probenames <- rownames(deTable[deTable$adj.P.Val <= pcut,])
  probelist <- str_split(probenames, pattern = ":")
  probenames <- unlist(lapply(probelist, function(a) a[5]))
  probes <- hm450[probenames]
  genes <- getNearestGene(probes)  
  
  # clean up #
  geneSymbols <- genes$nearestGeneSymbol
  geneSymbols[is.na(geneSymbols)] <- "none"
  geneSymbols[geneSymbols == ""] <- "none"
  
  # convert to entrez #
  geneEntrez <- mget(org.Hs.egALIAS2EG, x=geneSymbols, ifnotfound=NA) 
  geneEntrez <- unlist(lapply(geneEntrez, function(a) {
    if (length(a) > 1) {
      str_join(a, collapse=";")
    } else {
      if (is.na(a)) {
        "NA"
      } else {
        a
      }
    }
  }))
  tab1 <- cbind(genes, as.data.frame(seqnames(probes)), data.frame(Pos=probes$probeStart, Entrez=geneEntrez))
  tab1
}

mapToGenes2 <- function(deTable, pcut=0.05)
{
  require(stringr)  
  library("IlluminaHumanMethylation450k.db")
  probenames <- rownames(deTable[deTable$adj.P.Val <= pcut,])
  probelist <- str_split(probenames, pattern = ":")
  probenames <- unlist(lapply(probelist, function(a) a[5]))
  
  symbols  <- mget(probenames, IlluminaHumanMethylation450kSYMBOL, ifnotfound = NA)
  chrs     <- mget(probenames, IlluminaHumanMethylation450kCHR37, ifnotfound = NA)
  chrpos   <- mget(probenames, IlluminaHumanMethylation450kCHRLOC, ifnotfound = NA)
  probeloc <- mget(probenames, IlluminaHumanMethylation450kPROBELOCATION, ifnotfound = NA)
  coord    <- mget(probenames, IlluminaHumanMethylation450kCPGCOORDINATE, ifnotfound = NA)
  cgploc   <- mget(probenames, IlluminaHumanMethylation450kCPGILOCATION, ifnotfound=NA)
  data.frame(Genes=symbols, Chr=chrs, Pos=probeloc)
}


