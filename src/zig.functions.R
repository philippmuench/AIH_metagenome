# this function will create a fitZig model based on a biom table and a given sample grouping
getZigModBasedOnGrouping <- function(biom, grouping){
  mod <- model.matrix(~ 0 + grouping)
  settings <- zigControl(maxit = 10, verbose = TRUE)
  res <- fitZig(obj = biom, mod = mod, control = settings)
  return(res)
}

# fit model with contrast and return all significant ones (FDR corrected p-value < 0.05), add tax annotation
fitModelWithContrast <- function(biom, zigfit, contrast, onlysig=T, annotation=T, phylum =F){
	# fit model
	fit <- contrasts.fit(zig.treatment.fit , contrast)
	fit.ebayes <- eBayes(fit)
	# get best ones and make FDR correction
	top <- topTable(fit.ebayes, adjust.method="none", number = 100000)
	top$adj.P.Val <- p.adjust(top$adj.P.Val , method="fdr")
	if (onlysig){
		top <- top[which(top$adj.P.Val < 0.05),] # remove non-significant ones
	}
	# annotate taxa to top table
	if (annotation){
	 	if (phylum){
		  	# for the manual created MRexperiment we just have phylum lvl annotation
		 	otu.num <- match(rownames(top), rownames(biom))
			annot <- fData(biom)
			annot.phylum <- as.matrix(annot[match(rownames(top), rownames(annot)),]$V2)
			top <- cbind(top, annot.phylum, otu.num)
		} else {
			otu.num <- match(rownames(top), rownames(biom))
			annot <- fData(biom)
			annot.x2 <- as.matrix(annot[match(rownames(top), rownames(annot)),]$X2)
			annot.x3 <- as.matrix(annot[match(rownames(top), rownames(annot)),]$X3)
			annot.x4 <- as.matrix(annot[match(rownames(top), rownames(annot)),]$X4)
			annot.x5 <- as.matrix(annot[match(rownames(top), rownames(annot)),]$X5)
			annot.x6 <- as.matrix(annot[match(rownames(top), rownames(annot)),]$X6)
			annot.x7 <- as.matrix(annot[match(rownames(top), rownames(annot)),]$X7)
			top <- cbind(top,annot.x2, annot.x3, annot.x4, annot.x5, annot.x6, annot.x7, otu.num)
		}
	}
 	return(top)
}

# this function generates a volcano plot (logFC vs logP)
createVolcan <- function(top.contrast, avExprAsSize=F, phylum=F, title="Volcano plot"){
	require(ggplot2)
 	df <- as.data.frame(top.contrast)
 	if(phylum){
 		df$p <- df$annot.phylum
 	} else {
	 	split <- strsplit(as.character(as.matrix(df$annot.x2)),
						"_", fixed = TRUE)
	 	df$p <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 1)
	 	split <- strsplit(as.character(as.matrix(df$annot.x3)),
						"_", fixed = TRUE)
	 	df$c <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 1)
	 	split <- strsplit(as.character(as.matrix(df$annot.x4)),
						"_", fixed = TRUE)
		df$o <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 1)
		split <- strsplit(as.character(as.matrix(df$annot.x5)),
						"_", fixed = TRUE)
		df$f <- sapply(strsplit(sapply(split, "[", 3),"_"),"[", 1)
		df <- df[-which(is.na(df$p)),] # remove all without phylum annotation
		df <- df[-which(is.na(df$o)),] # remove all without order annotation
		df$col <- taxCol(df$p, df$c) # assign colors to data
 	}
	df$adj.P.Val <- log10(df$adj.P.Val)
	df.annot <- df[which(df$adj.P.Val < -3),]
	cat(paste(nrow(df.annot), "are significant"))
	if (avExprAsSize){
		c <- ggplot(df, aes(x=logFC, y=adj.P.Val, size=AveExpr))
	} else {
		c <- ggplot(df, aes(x=logFC, y=adj.P.Val, size=0.5))
	}

	

# generates for on otu a plot
plotOTU2 <- function (obj, otu, classIndex, log = TRUE, norm = TRUE, jitter.factor = 1,
											pch = 21, labs = TRUE, xlab = NULL, ylab = NULL, jitter = TRUE,
											...){

	mat = returnAppropriateObj(obj, norm, log)
	l = lapply(classIndex, function(j) {
		mat[otu, j]
	})
	z = posteriorProbs(obj)
	y = unlist(l)
	x = rep(seq(along = l), sapply(l, length))
	if (!is.null(z)) {
		z = 1 - z
		lz = lapply(classIndex, function(j) {
			(z[otu, j])
		})
		z = unlist(lz)
		blackCol = t(col2rgb("black"))
		col = rgb(blackCol, alpha = z)
	}
	else {
		blackCol = t(col2rgb("black"))
		col = rgb(blackCol)
	}
	if (jitter)
		x = jitter(x, jitter.factor)
	if (is.null(ylab)) {
		ylab = "Normalized log(cpt)"
	}
	if (is.null(xlab)) {
		xlab = "Groups of comparison"
	}

	#print(  as.matrix(x))
	#print(  as.matrix(x))
	#  p <- ggplot(ToothGrowth, aes(x=dose, y=len)) +
	#    geom_boxplot()

	plot(x, y, col = col, pch = pch, bg = col, xlab = xlab, ylab = ylab,
			 xaxt = "n", ...)
	if (labs == TRUE) {
		gp = names(classIndex)
		axis(1, at = seq(1:length(gp)), gp)
	}
	invisible(list(x = x, y = y))
}

# generates OTU plots for a annotated top file generated with fitModelWithContrast()
plotMultipleOtus <- function(biom.table, top, outputfile="figures/otu_plot.pdf"){
	require(metagenomeSeq)
	# iterate over lsit of OTUs and generate a OTU plot for each one and save it in one pdf file
	pdf(outputfile, height = 5, width = 10)
	for (otu in top$otu.num){
		# get the annotation of this otu
		otu.tax <- paste(otu,top[match(otu, top$otu.num),]$annot.x4,
										 top[match(otu, top$otu.num),]$annot.x5,
										 top[match(otu, top$otu.num),]$annot.x6,
										 top[match(otu, top$otu.num),]$annot.x7
		)
		# generate plot with annotation as title
		plotOTU2(biom.table, otu = otu, classIndex, main = otu.tax, jitter.factor = 1, jitter=T)
	}
	dev.off()
}
	c <- c + geom_point(alpha=1) + guides(size = FALSE)
	c <- c + scale_y_reverse()
	c <- c + ggtitle(title)
	#c <- c + scale_color_manual(values = df$col)
	c <- c +geom_hline(yintercept = -1.301) + theme_minimal() # + guides(color = FALSE)
	if(phylum){
		c <- c + geom_text_repel(data=df, aes(x=logFC,y=adj.P.Val, label=p, size=0.1))
	} else {
		c <- c + geom_text_repel(data=df.annot, aes(x=logFC,y=adj.P.Val, label=f, size=1))
	}
return(c)
}
