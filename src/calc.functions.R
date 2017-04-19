#' by Philipp C. Münch (philipp.muench@helmholtz-hzi.de)
#' this function will remove a Treatment type from the biom table (requires that there is a Treatment column)
removeSamplefromBiom <- function(biom.table, sample){
  require(metagenomeSeq)
  cat("remove sample")
  biom.table.p <- pData(biom.table) # get pdata annotation
  samplesToRemove <- which(biom.table.p$Treatment == sample)
  biom.table.reduced <- biom.table[, -samplesToRemove] # generate new obj
  return(biom.table.reduced)
}

# this function will remove features from a biom table that are present in < threshold samples
removeFeaturesbasedOnMRcounts <- function(biom.table, threshold){
  require(metagenomeSeq)
  cat("remove rare features")
  rareFeatures <- which(rowSums(MRcounts(biom.table) > 0) < threshold)
  biom.table.reduced <- biom.table[-rareFeatures, ]
  return(biom.table.reduced)
}

# this function will normalize a biom table
normalizeBiomTable <- function(biom.table){
  require(metagenomeSeq)
  cat("normalize OTU table")
  p <- cumNormStat(biom.table, pFlag = TRUE, main = " OTU data")
  biom.table.normalized <- cumNorm(biom.table, p = p)
  return(biom.table.normalized)
}
removeRelAb <- function(data.mat,
                        threshold = 0.01){ # remove genera with leth than (1%, threshold% as their maximum relative abundance
  # determine the maximum relative abundance for each column
  maxab <- apply(data.mat, 2, max)
  # find out which genera to remove
  n1 <- names(which(maxab < threshold))
  cat(paste("removed", length(n1),"taxa (lower than", threshold*100,"%)" , "\n"))
  # new matrix without genera
  data.mat.filtered <- data.mat[, -which(colnames(data.mat) %in% n1)]
  return(data.mat.filtered)
}

makeOtu <- function(raw.data, method="proportion"){
  # basic manipulations on raw OTU data
  rownames(raw.data) <- as.character(as.matrix(raw.data$ID)) #add names to row
  out.data <- raw.data[,-1] #remove id lable
  out.data$other <- NULL # remove taxa "other"
  out.data <- t(out.data) #swap row/col
  if (method == "proportion") {
    #transform the raw counts of reads to proportions within a sample
    cat(paste("transform the raw counts to proportions", "\n"))
    out.data <- out.data/rowSums(out.data)
    # todo: FC transform
    # todo(pmuench): log(FC) transform
    #data.prop.1 <- data.prop.1 - data.prop.1[which(rownames(data.prop.1) == "control"),]
  }
  return(out.data)
}

processTaxaID <- function(t.data){
  # process the ID tax information
  t.data.tax <- as.data.frame(t.data)
  split <- strsplit(as.character(as.matrix(rownames(t.data.tax))),
                    ";", fixed = TRUE)
  if (sapply(split, length)[1] == 2){ # L2 otu table
    cat("process taxa (L2)")
    t.data.tax$k <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
    t.data.tax$p <- sapply(strsplit(sapply(split, "[", 2),"_"),"[", 3)
  } else if (sapply(split, length)[1] == 3){ # L3 otu table
    cat("process taxa (L3)")
    t.data.tax$k <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
    t.data.tax$p <- sapply(strsplit(sapply(split, "[", 2),"_"),"[", 3)
    t.data.tax$c <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 3)
  } else if (sapply(split, length)[1] == 4){ # L4 otu table
    cat("process taxa (L4)")
    t.data.tax$k <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
    t.data.tax$p <- sapply(strsplit(sapply(split, "[", 2),"_"),"[", 3)
    t.data.tax$c <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 3)
    t.data.tax$o <- sapply(strsplit(sapply(split, "[", 4),"_"),"[", 3)
  } else if (sapply(split, length)[1] == 5){ # L5 otu table
    cat("process taxa (L5)")
    t.data.tax$k <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
    t.data.tax$p <- sapply(strsplit(sapply(split, "[", 2),"_"),"[", 3)
    t.data.tax$c <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 3)
    t.data.tax$o <- sapply(strsplit(sapply(split, "[", 4),"_"),"[", 3)
    t.data.tax$f <- sapply(strsplit(sapply(split, "[", 5),"_"),"[", 3)
  } else if (sapply(split, length)[1] == 7){ # L7 otu table
    cat("process taxa (L7)")
    t.data.tax$k <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
    t.data.tax$p <- sapply(strsplit(sapply(split, "[", 2),"_"),"[", 3)
    t.data.tax$c <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 3)
    t.data.tax$o <- sapply(strsplit(sapply(split, "[", 4),"_"),"[", 3)
    t.data.tax$f <- sapply(strsplit(sapply(split, "[", 5),"_"),"[", 3)
    t.data.tax$g <- sapply(strsplit(sapply(split, "[", 6),"_"),"[", 3)
    t.data.tax$s <- sapply(strsplit(sapply(split, "[", 7),"_"),"[", 3)
  } else if (sapply(split, length)[1] == 6){ # L6 otu table
    cat("process taxa (L6)")
    t.data.tax$k <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
    t.data.tax$p <- sapply(strsplit(sapply(split, "[", 2),"_"),"[", 3)
    t.data.tax$c <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 3)
    t.data.tax$o <- sapply(strsplit(sapply(split, "[", 4),"_"),"[", 3)
    t.data.tax$f <- sapply(strsplit(sapply(split, "[", 5),"_"),"[", 3)
    t.data.tax$g <- sapply(strsplit(sapply(split, "[", 6),"_"),"[", 3)
  } else {
    cat("error")
  }
  return(t.data.tax)
}

processTaxaDiff <- function(otu.data,
                            column.control = "control",
                            column.normalize = c("non", "AIH")){
  # caluclate difference to sample in % of abundance
  column.control.id <- which(colnames(otu.data) == column.control)

  for (column.normalize.element in column.normalize){
    column.normalize.id <- which(colnames(otu.data) == column.normalize.element)
    otu.data[,column.normalize.id] <- otu.data[,column.normalize.id] - otu.data[,column.control.id] # susbtract control from list of column.normalize
  }
  otu.data[,column.control.id] <- otu.data[,column.control.id] - otu.data[,column.control.id] # shoud be 0
  return(otu.data)
}

# this function will remove a sample from a distance file
removeSampleFromDist <- function(dist, sampleName="HIV", mapping){
  idToRemove <- which(mapping == sampleName)
  dist.clean <- dist[-idToRemove, -idToRemove]
  cat(paste("removed",length(idToRemove), "samples from distance matrix with sample name", sampleName))
  return(dist.clean)
}
