getOtuCorPvalue <- function(otu, comparison){
corMat <- data.frame(name=colnames(otu), cor=rep(0, ncol(otu)), pval=rep(-1, ncol(otu)))
  for (i in 1:ncol(otu)){
    # make get cor pvalues
    corMat[i,]$pval <- cor.test(otu[,i],comparison, alternative="two.sided")$p.value
    corMat[i,]$cor <- cor.test(otu[,i],comparison, alternative="two.sided")$estimate
  }
  corMat$p.adj <- p.adjust(corMat$pval, method="fdr")
  return(corMat)
}

drawCorMarker <- function(df){
  require(ggplot2)
  a <- ggplot(df, aes(reorder(annot.name,p.adj), cor, color=type))
  a <- a + stat_summary(geom = 'errorbar', fun.data = 'seFunc', width = 0, show_guide = F, position="dodge")
  a <- a + stat_summary(geom = 'point', fun.y = 'mean', size = 2, shape = 21,position="dodge") + theme_minimal()
  a <- a + theme(axis.text.x=element_text(angle = -45, hjust = 0))
  a <- a + geom_text(aes(label=is.sig), size=8)
  a <- a + scale_color_manual(values = c("black","grey50", "#00c094ff"),labels = c("AIH","control", "healthy"))
  return(a)  
}

processTaxa <- function(df){
  split <- strsplit(as.character(as.matrix(df$species)),
                      ";", fixed = TRUE)
  split2 <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
  df$species.split <- split2
  split <- strsplit(as.character(as.matrix(df$genus)),
                      ";", fixed = TRUE)
  split2 <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
  df$genus.split <- split2
  df$genus.first <- substr(df$genus.split, 1,1)
  df$annot.name <- paste(df$genus.first, ". ", df$species.split, sep="")
  return(df)
}

drawCorMarkerOtu <- function(df, title="correlation markers to alpha diversity", all=T){
	a <- ggplot(df, aes(x=adiv, y=value)) 
	a <- a + geom_smooth(method=lm,se=FALSE, aes(color=type)) 
	if (all){
		a <- a + geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black")		
	}
	a <- a + geom_point(shape=1, aes(color=type))
	a <- a + facet_grid(variable ~ ., scales="free_y") + theme_bw()
	a <- a + scale_color_manual(values = c("black","#00c094ff", "grey50", "darkblue"),labels = c("AIH","healthy", "control", "all"))
	a <- a + ggtitle(title) 
	return(a)
}

printCorrelationCoef <- function(df.adiv){
	# control
	cat("control\n")
	print(cor.test(df.adiv[which(df.adiv$type =="others"),]$adiv, df.adiv[which(df.adiv$type =="others"),]$ifap, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="others"),]$adiv, df.adiv[which(df.adiv$type =="others"),]$lbp, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="others"),]$adiv, df.adiv[which(df.adiv$type =="others"),]$sCD14, alternative="two.sided"))
	# healthy
	cat("healthy\n")
	print(cor.test(df.adiv[which(df.adiv$type =="control"),]$adiv, df.adiv[which(df.adiv$type =="control"),]$ifap, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="control"),]$adiv, df.adiv[which(df.adiv$type =="control"),]$lbp, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="control"),]$adiv, df.adiv[which(df.adiv$type =="control"),]$sCD14, alternative="two.sided"))
	# AIH
	cat("AIH\n")
	print(cor.test(df.adiv[which(df.adiv$type =="AIH"),]$adiv, df.adiv[which(df.adiv$type =="AIH"),]$ifap, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="AIH"),]$adiv, df.adiv[which(df.adiv$type =="AIH"),]$lbp, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="AIH"),]$adiv, df.adiv[which(df.adiv$type =="AIH"),]$sCD14, alternative="two.sided"))

	cat("all\n")
	print(cor.test(df.adiv[which(df.adiv$type =="AIH" | df.adiv$type =="control" | df.adiv$type =="others" ),]$adiv, df.adiv[which(df.adiv$type =="AIH" | df.adiv$type =="control" | df.adiv$type =="others"),]$ifap, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="AIH" | df.adiv$type =="control" | df.adiv$type =="others"),]$adiv, df.adiv[which(df.adiv$type =="AIH" | df.adiv$type =="control" | df.adiv$type =="others"),]$lbp, alternative="two.sided"))
	print(cor.test(df.adiv[which(df.adiv$type =="AIH" | df.adiv$type =="control" | df.adiv$type =="others"),]$adiv, df.adiv[which(df.adiv$type =="AIH" | df.adiv$type =="control" | df.adiv$type =="others"),]$sCD14, alternative="two.sided"))
}

matchMarker <- function(adiv, marker){
	adiv$ifap <- marker[match(adiv$sample, marker$id),]$ifabp
	adiv$lbp <- marker[match(adiv$sample, marker$id),]$lbp
	adiv$sCD14 <- marker[match(adiv$sample, marker$id),]$sCD14
	adiv$adiv <- adiv$value
	adiv$value <- NULL
	return(adiv)
}
