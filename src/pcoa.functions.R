# code adapted from http://enterotype.embl.de/enterotypes.html
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

# draws pcoa using ggbiplot function
drawPcoaGbiplot <- function(fit, grouping){
  require(ggbiplot)
  c <- ggbiplot(fit,obs.scale = 1, var.axes=F,var.scale=1,
                groups=grouping,ellipse=T,circle=T,varname.size=3)
  c <- c + theme_bw()
  c <- c + scale_color_manual(values=c("darkred", "darkgreen", "darkblue"))
  return(c)
}

# calculates the explained prob. for pcoa plot
pcoaProp <- function(obj){
  eig_tot <- sum(obj$values$Eigenvalues)
  prob <-  obj$values$Eigenvalues / eig_tot
  return(prob)
}

# draw pcoa using ggplot
drawPcoa <- function(obs.pcoa, df, sample.names=sample.type, colorvalues=colorvalues, label=F)
{
  p <- ggplot(df, aes(x=MDS1,y=MDS2))
  p <- p + geom_point(aes(size=4, alpha=0.8, colour= sample.names, label =sample.names)) # draw points
  #  if (ellipse){
  #   p <- p + geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y), color="grey80")
  #  p <- p + geom_path(data=conf.rgn) + theme_bw() + geom_hline(yintercept=0,linetype=3)
  # p <- p + geom_vline(xintercept=0,linetype=2) + ggtitle("")
  # p <- p + geom_point(data=centroids, size=1, color="grey80") # draw centroids
  # }
  
  p <- p + xlab(paste("PCo 1 (", format(pcoaProp(obs.pcoa)[1]*100, digits=4), "%)",  sep=""))
  p <- p + ylab(paste("PCo 2 (", format(pcoaProp(obs.pcoa)[2]*100, digits=4),"%)", sep=""))
  p <- p + theme(legend.position="none")
  p <- p + theme(plot.title = element_text(size = 8),
                 axis.title.x = element_text(size = 8),
                 axis.title.y = element_text(size = 8))
  p <- p + theme(plot.title = element_text(size = 8),
                 axis.title.x = element_text(size = 8),
                 axis.title.y = element_text(size = 8))
  p <- p + scale_color_manual(values = colorvalues)
  p <- p + theme_bw() + geom_hline(yintercept = 0,linetype = 3)
  p <- p + geom_vline(xintercept = 0,linetype = 3)
  if (label){
    p <- p + geom_label(size = 2, aes(fill = type.names), colour = "white", fontface = "bold")
  }
  return(p)
}


# rubens unconstrained pcoa function
plot_pcoa <- function(p, d, variables, col_var, pch_var, col_comp, pch_comp,
                      shapes, colors, file_name, svg=F, pcoa_var, cex.main=0.8, cex=0.8) {
  
  col <- rep("black", dim(p)[1])
  pch <- rep(19, dim(p)[1])
  for (i in 1:length(col_var)){
    index <- grep(col_var[i], d[rownames(p), col_comp[i]])
    col[index] <- colors[i]
  }
  for (i in 1:length(pch_var)){
    index <- grep(pch_var[i], d[rownames(p), pch_comp[i]])
    pch[index] <- shapes[i]
  }
  xlab <- paste("PCo 1 (", format(pcoa_var[1]*100, digits=4), "%)", sep="")
  ylab <- paste("PCo 1 (", format(pcoa_var[2]*100, digits=4), "%)", sep="")
  main <- "Unconstrained PCoA"
  if(svg) svg(file=file_name) else pdf(file=file_name)
  plot(p, col=col, pch=pch, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex=cex)
  abline(v=0, h=0, lty="dotted")
  dev.off()
  
}