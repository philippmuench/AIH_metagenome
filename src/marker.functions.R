#' by Philipp C. Münch (philipp.muench@helmholtz-hzi.de)

getOtuCorPvalueMulti <- function(otu, comparisons, pval=FALSE, totalcomparisons){
 # corMat <- data.frame(name=colnames(otu))
corMat <-   matrix(, nrow = ncol(otu), ncol =ncol(comparisons))
colnames(corMat) <- colnames(comparisons)
rownames(corMat) <- colnames(otu)
  for (column in 1:ncol(comparisons)){

    cat(paste(colnames(comparisons)[column], "\n"))
     for (i in 1:ncol(otu)){
      if (pval){
    corMat[i,column] <- p.adjust(cor.test(otu[,i],as.numeric(as.matrix(comparisons[,column])), alternative="two.sided")$p.value, n=totalcomparisons, method="fdr")
      } else {
           corMat[i,column] <- cor.test(otu[,i],as.numeric(as.matrix(comparisons[,column])), alternative="two.sided")$estimate
      }

    }
  }
  return(corMat)
}
