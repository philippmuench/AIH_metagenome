uniqueFeatures<-function(obj,cl,nsamples=0,nreads=0){
  if (class(obj) == "MRexperiment") {
    mat = MRcounts(obj, norm = FALSE, log = FALSE)
  }
  else if (class(obj) == "matrix") {
    mat = obj
  }
  else {
    stop("Object needs to be either a MRexperiment object or matrix")
  }
  res = by(t(mat),cl,colSums)
  res = do.call("rbind",res)
  kreads = (colSums(res==0)>0)

  mat = mat>0
  resPos = by(t(mat),cl,colSums)
  resPos = do.call("rbind",resPos)
  ksamples = (colSums(resPos==0)>0)

  featureIndices = intersect(which(ksamples),which(kreads))
  numberReads = t(res[,featureIndices])
  colnames(numberReads) = paste("Reads in",colnames(numberReads))
  numberPosSamples = t(resPos[,featureIndices])
  colnames(numberPosSamples) = paste("Samp. in",colnames(numberPosSamples))
  featureIndices = featureIndices
  featureNames = rownames(mat[featureIndices,])





  df = cbind(featureIndices,numberPosSamples,numberReads)
  interesting = which(rowSums(numberReads)>=nreads & rowSums(numberPosSamples)>=nsamples)
  df[interesting,]
}
