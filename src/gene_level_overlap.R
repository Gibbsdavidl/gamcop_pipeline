source("src/gene_list_overlap.R")
add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}

subset_de_table <- function(deTable, columnsToInclude){
  featureNames <- rownames(deTable)
  tab <- deTable[, columnsToInclude]
  tab <- cbind(featureNames, tab)
  return(tab)
}

# output overlap statistics
# gene1 and gene 2 - character strings that represent the gene identifier to map between tables
# comparsion e.g. "expression_methylation"
generate_de_overlap_table <- function(deTable1, deTable2){
  col1 <- grep("gene", names(deTable1), ignore.case = T)
  col2 <- grep("gene", names(deTable2), ignore.case = T)
  table <- merge(deTable1, deTable2, by.x = col1, by.y=col2)
  return(table)
}

make_overlap_tables <- function(mRNA, miRNA, Meth, targetPhenotypeName, writingDir){
  sigColumns <- c("logFC", "adj.P.Val")
  
  columnsToInclude_mRNA <- c(sigColumns, "GeneName")
  rna <- subset_de_table(mRNA, columnsToInclude_mRNA)
  
  columnsToInclude_microRNA <-c(sigColumns, "HGNC.gene")
  mirna <- subset_de_table(miRNA, columnsToInclude_microRNA)
  
  columnsToInclude_meth <- c(sigColumns, "nearestGeneSymbol")
  meth <-  subset_de_table(Meth, columnsToInclude_meth)
  
  if(!is.null(mRNA) & !is.null(miRNA)){
    root <-  add_date_tag(paste("combined_DE_table", targetPhenotypeName, c("mRNA_miRNA"), sep = "_"), ".txt")
    filename_path <- paste0(writingDir, root)
    res <- generate_de_overlap_table(rna, mirna)
    names(res) <- gsub(pattern = ".x",replacement = "_mRNA", x = names(res))
    names(res) <- gsub(pattern = ".y",replacement = "_miRNA", x = names(res))    
    write.table(res, file = filename_path, sep = "\t", row.names = F, quote = F)
  }  
  
  if(!is.null(mRNA) & !is.null(Meth)){
    root <-  add_date_tag(paste("combined_DE_table", targetPhenotypeName, c("mRNA_meth"), sep = "_"), ".txt")
    filename_path <- paste0(writingDir, root)
    res <- generate_de_overlap_table(rna, meth)
    names(res) <- gsub(pattern = ".x",replacement = "_mRNA", x = names(res))
    names(res) <- gsub(pattern = ".y",replacement = "_meth", x = names(res))    
    write.table(res, file = filename_path, sep = "\t", row.names = F, quote = F)
  }
  
  
  if(!is.null(miRNA) & !is.null(Meth)){
    root <-  add_date_tag(paste("combined_DE_table", targetPhenotypeName, c("miRNA_meth"), sep = "_"), ".txt")
    filename_path <- paste0(writingDir, root)
    res <- generate_de_overlap_table(mirna, meth)
    names(res) <- gsub(pattern = ".x",replacement = "_miRNA", x = names(res))
    names(res) <- gsub(pattern = ".y",replacement = "_meth", x = names(res))
    write.table(res, file = filename_path, sep = "\t", row.names = F, quote = F)
  }  
}