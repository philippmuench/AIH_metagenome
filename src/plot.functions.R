#' by Philipp C. Münch (philipp.muench@helmholtz-hzi.de)

require(ggplot2)

# make a standard error function for plotting
seFunc <- function(x){
 se <- sd(x) / sqrt(sum(!is.na(x)))
 lims <- c(mean(x) + se, mean(x) - se)
 names(lims) <- c('ymin', 'ymax')
 return(lims)
}


mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line


# this function generates color pattern for taxonomic groups
taxCol <- function(lvl1, lvl2, rand=F){
	require(RColorBrewer)
	obj <- data.frame(lvl1, lvl2)
	nr <- nrow(obj)
	tbl <- with(obj, table(lvl1)[order(unique(lvl1))])
	# generate a color lookup table
	pNames <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Tenericutes",  "Verrucomicrobia", "Fusobacteria", "Cyanobacteria", "Deferribacteres")
	pColor <-  c("darkorange4", "aquamarine4", "chartreuse4", "bisque4", "blue4", "brown4", "yellow", "orange", "brown")
	colorLookup <- data.frame(name=pNames, color=pColor)
	cols <- colorLookup[match(names(tbl), colorLookup$name),]$color
	# get base color from lookup table
	if (rand){
		# assign colors in a random manner
		qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] # get quantitaive colors
		col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
		cols <- sample(col_vector, length(unique(lvl1)))
	}
	#cols <- c('cyan2','red','orange','green','dodgerblue2', 'black', 'yellow') # assign here a list of colors
	cols <- unlist(Map(rep, cols, tbl))
	col1 <- rep(NA, nr)
	col2 <- rep(NA, nr)
	for (i in 1:nr) {
		## create color/shades
		rgb <- col2rgb(cols[i])
		col2[i] <- rgb(rgb[1], rgb[2], rgb[3], 190 / sequence(tbl)[i], maxColorValue = 255)
		## repeat above for the main groups
		rgb <- col2rgb(cols[i])
		col1[i] <- rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
	}
	return(col2)
}

# this function generates color pattern based on phylum assignment
taxColPhylum <- function(lvl1){
	require(RColorBrewer)
	obj <- data.frame(lvl1)
	nr <- nrow(obj)
	tbl <- with(obj, table(lvl1)[order(unique(lvl1))])
	# generate a color lookup table
	pNames <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Tenericutes",  "Verrucomicrobia", "Fusobacteria", "Cyanobacteria", "Deferribacteres")
	pColor <-  c("darkorange4", "aquamarine4", "chartreuse4", "bisque4", "blue4", "brown4", "yellow", "orange", "brown")
	colorLookup <- data.frame(name=pNames, color=pColor)
	cols <- colorLookup[match(names(tbl), colorLookup$name),]$color
	return(cols)
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



# draw tri plot with log normalized otu information and mark significant ones based on two input tables "sigTable1" and "sigTable2"
triPlotOTU <- function(biom.table, sigTable1, sigTable2, table1_name, table2_name){
	require(metagenomeSeq)
	require(ggplot2)
	require(ggtern)
	mat <- returnAppropriateObj(biom.table, norm=T, log=T)
	type <- pData(biom.table)$Treatment
	type[which(type == "control")]  <- "healthy"
	type[which(type != "AIH" & type != "healthy")] <- "control"
	type2 <- rep("NA", nrow(mat))
	type2[which(!is.na(match(rownames(mat), rownames(sigTable1))))] <- table1_name
	type2[which(!is.na(match(rownames(mat), rownames(sigTable2))))] <- table2_name
	mat.healthy <- rowMeans(mat[,which(type == "healthy")])
	mat.AIH <- rowMeans(mat[,which(type == "AIH")])
	mat.control <- rowMeans(mat[,which(type == "control")])
	df <- as.data.frame(cbind(mat.healthy, mat.AIH, mat.control))
	df$average <- rowMeans(df)
	df$type2 <- type2
	# add taxa to significant ones
	df$x2 <- "NA"
	df$x2[match(rownames(sigTable1), rownames(df))] <- as.character(sigTable1$annot.x2)
	df$x3 <- "NA"
	df$x3[match(rownames(sigTable1), rownames(df))] <- as.character(sigTable1$annot.x3)
	df$x4 <- "NA"
	df$x4[match(rownames(sigTable1), rownames(df))] <- as.character(sigTable1$annot.x4)
	df$x5 <- "NA"
	df$x5[match(rownames(sigTable1), rownames(df))] <- as.character(sigTable1$annot.x5)
	df$pval <- "NA"
	df$pval[match(rownames(sigTable1), rownames(df))] <- -log10(as.numeric(sigTable1$adj.P.Val))
	df$pval[match(rownames(sigTable2), rownames(df))] <- -log10(as.numeric(sigTable2$adj.P.Val))
	df.clean <- df[-which(df$x3=="NA"),] # remove all without tax annotation
	# plot
	g <- ggtern(df.clean,aes(mat.healthy,mat.AIH, mat.control, color=x3))
	#g <- g + geom_point(aes(size = as.numeric(healthy.pval)), alpha=0.5) # add point AIH
	g <- g + geom_point(aes(size = average), alpha=0.5) # add point AIH
	g <- g + facet_grid(.~ type2)
	g <- g + theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="#FFFFFF"))
	g <- g +theme_minimal()# + geom_text(aes(label=x3))
	return(g)
}

# draw list of taxa as bar chart
pieChartTaxaList <- function(listOfTaxa){
	require(ggplot2)
	occ <- table(listOfTaxa)
	occ.m <- melt(occ)
	g <- ggplot(occ.m,aes("",value, fill=listOfTaxa))
	g <- g + geom_bar(stat = "identity") + coord_polar("y")
	return(g)
}

givemedonutsorgivemedeath <- function(file, width = 15, height = 11) {
	## house keeping
	if (missing(file)) file <- getwd()
	plot.new(); op <- par(no.readonly = TRUE); on.exit(par(op))

	pdf(file, width = width, height = height, bg = 'snow')

	## useful values and colors to work with
	## each group will have a specific color
	## each subgroup will have a specific shade of that color
	nr <- nrow(browsers)
	width <- max(sqrt(browsers$share)) / 0.8

	tbl <- with(browsers, table(browser)[order(unique(browser))])
	cols <- c('cyan2','red','orange','green','dodgerblue2', 'black', 'yellow')
	cols <- unlist(Map(rep, cols, tbl))

	## loop creates pie slices
	plot.new()
	par(omi = c(0.5,0.5,0.75,0.5), mai = c(0.1,0.1,0.1,0.1), las = 1)
	for (i in 1:nr) {
		par(new = TRUE)

		## create color/shades
		rgb <- col2rgb(cols[i])
		f0 <- rep(NA, nr)
		f0[i] <- rgb(rgb[1], rgb[2], rgb[3], 190 / sequence(tbl)[i], maxColorValue = 255)

		## stick labels on the outermost section
		lab <- with(browsers, sprintf('%s: %s', version, share))
		if (with(browsers, share[i] == max(share))) {
			lab0 <- lab
		} else lab0 <- NA

		## plot the outside pie and shades of subgroups
		pie(browsers$share, border = NA, radius = 5 / width, col = f0,
				labels = lab0, cex = 1.8)

		## repeat above for the main groups
		par(new = TRUE)
		rgb <- col2rgb(cols[i])
		f0[i] <- rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)

		pie(browsers$share, border = NA, radius = 4 / width, col = f0, labels = NA)
	}

	## extra labels on graph

	## center labels, guess and check?
	text(x = c(-.05, -.05, 0.15, .25, .3), y = c(.08, -.12, -.15, -.08, -.02),
			 labels = unique(browsers$browser), col = 'white', cex = 1.2)
	dev.off()
}

assignColors <- function(lower.level # lower level tax. assignment e.g. phyla
												 , higher.level){ # higher level tax. assignment e.g. class
	# assign colors according to taxonomic assignment, same taxa on lower.level levl
	# will have similar color, and taxa on different lower.levels have different colors
	num.col.groups <- length(unique(lower.level)) # num of different color groups
	colors <- rep("#FFFFFF", length(higher.level))
	hue_end <- 360*(num.col.groups-1)/num.col.groups # get max. hue
	temp_start <- 0
	ind <- 1
	for (tax in factor(lower.level)){
		alpha.values <- 1 - seq(from = 0, to = 0.4,
														length.out = length(which(lower.level == tax)))
#		increment <- hue_end/num.col.groups # possible color range for one group
#		colors[which(lower.level == tax)] <- rainbow_hcl(length(which(lower.level == tax)),)temp_start <-  temp_start + increment # update temp start
#		ind <- ind + 1
	}
	return(colors)
}

# generates are barplot for a given data frame
barplotClass <- function(df, file="figures/barplot_class.pdf", tooSmall = 0.2){
	require(ggplot2)
	print(paste("plot written to", file))
	df$mean <- rowMeans(cbind(df$healthy,df$control,df$AIH)) # mean values over samples
	df$type <- "NA"
	df[which(df$mean < tooSmall),]$type <- "(a) low abundant"
	df[which(df$mean >= tooSmall),]$type <- "(b) high abundant"
	df$mean <- NULL
	df.m <- melt(df)
	c <- ggplot(df.m, aes(x=variable, y=value, fill=paste(p,c)))
	c <- c + geom_bar(stat="identity", width = 0.4)
	c <- c + theme(legend.position="right")
	c <- c + scale_fill_manual(values = df.m$col)
	c <- c + theme_minimal() + theme(panel.grid.major = element_blank(),
																	panel.grid.minor = element_blank(),
																	panel.border = element_blank(),
																	panel.background = element_blank())
	c <- c + facet_grid(type~., scales="free", space="free_x", drop=T) #+ theme(panel.margin = unit(-1, "lines"))
	c <- c + xlab("") + ylab("")
	pdf(file, width = 5, height = 5)
	print(c)
	dev.off()
	return(c)
}

# generates are barplot for a given data frame
barplotPhylum<- function(df, file="figures/barplot_phylum.pdf", tooSmall = 0.2){
	require(ggplot2)
	print(paste("plot written to", file))
	df$mean <- rowMeans(cbind(df$healthy,df$control,df$AIH)) # mean values over samples
	df$type <- "NA"
	df[which(df$mean < tooSmall),]$type <- "(a) low abundant"
	df[which(df$mean >= tooSmall),]$type <- "(b) high abundant"
	df$mean <- NULL
	df.m <- melt(df)
	c <- ggplot(df.m, aes(x=variable, y=value, fill=p))
	c <- c + geom_bar(stat="identity", width = 0.4)
	c <- c + theme(legend.position="right")
	c <- c + scale_fill_manual(values = df.m$col)
	c <- c + theme_minimal() + theme(panel.grid.major = element_blank(),
	                                 panel.grid.minor = element_blank(),
	                                 panel.border = element_blank(),
						panel.background = element_blank())
	c <- c + facet_grid(type~., scales="free", space="free_x", drop=T) + theme(panel.margin = unit(-1, "lines"))
	c <- c + xlab("") + ylab("")
	pdf(file, width = 5, height = 5)
	print(c)
	dev.off()
	return(c)
}


barplotGenus <- function(df, file="figures/barplot_class.pdf", tooSmall = 0.2){
	require(ggplot2)
	print(paste("plot written to", file))
	df$mean <- rowMeans(cbind(df$healthy,df$control,df$AIH)) # mean values over samples
	df$type <- "NA"
	df[which(df$mean < tooSmall),]$type <- "(a) low abundant"
	df[which(df$mean >= tooSmall),]$type <- "(b) high abundant"
	df$mean <- NULL
	df.m <- melt(df)
	c <- ggplot(df.m, aes(x=variable, y=value, fill=paste(c,g)))
	c <- c + geom_bar(stat="identity", width = 0.4)
	c <- c + theme(legend.position="right")
	c <- c + scale_fill_manual(values = df.m$col)
	c <- c + theme_minimal() + theme(panel.grid.major = element_blank(),
																	panel.grid.minor = element_blank(),
																	panel.border = element_blank(),
																	panel.background = element_blank())
	c <- c + facet_grid(type~., scales="free", space="free_x", drop=T) + theme(panel.margin = unit(-1, "lines"))
	c <- c + xlab("") + ylab("")

	pdf(file, width = 15, height = 5)
	print(c)
	dev.off()
	return(c)
}

barplotOrder <- function(df.m, file="figures/barplot_order.pdf"){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 5, height = 5)
	c <- ggplot(df.m, aes(x=variable, y=value, fill=paste(p,o)))
	c <- c + geom_bar(stat="identity", width = 0.4)
	c <- c + theme(legend.position="right")
	c <- c + scale_fill_manual(values = df.m$col) # + coord_polar("y")
	c <- c + theme_minimal() + theme(panel.grid.major = element_blank(),
																	 panel.grid.minor = element_blank(),
																	 panel.border = element_blank(),
																	 panel.background = element_blank())
	print(c)
	dev.off()
	return(c)
}

plotDifferenceClass <- function(df, file="figures/difference_class.pdf", jColors){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	# add type to colors
	jColors2 <- c("black", "grey")
	names(jColors2) <- c("AIH", "control")
	jColors <- c(jColors, jColors2)
	c <- ggplot(df, aes(x=reorder(c, -value), y=value*100))
	c <- c + geom_line(aes(group=c, color=c), size=2)
	c <- c + geom_point(size=5, shape=73,aes(color=variable)) + coord_flip()
	c <- c + scale_color_manual(values = jColors)
#  c <- c + scale_color_manual(values = c(control="grey", AIH="black"))
	c <- c + geom_hline(yintercept = 0, linetype = 3, color="grey80")
	c <- c + theme(legend.position="right")
#  c <- c + scale_color_manual(values = c("darkgrey","black"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())
	c <- c  + ylab("difference in RA compared to healthy samples") +
		xlab("")
	print(c)
	dev.off()
	return(c)
}

plotDifferencePhylum <- function(df, jColors){
	require(ggplot2)
#	print(paste("plot written to", file))

	# add type to colors
	jColors2 <- c("black", "grey")
	names(jColors2) <- c("AIH", "control")
	jColors <- c(jColors, jColors2)
	c <- ggplot(df, aes(x=reorder(p, -value), y=value*100))
	c <- c + geom_line(aes(group=p, color=p), size=2)
	c <- c + geom_point(size=5, shape=73,aes(color=variable)) + coord_flip()
	c <- c + scale_color_manual(values = jColors)
  #    c <- c + scale_color_manual(values = c(control="grey", AIH="black"))
	c <- c + geom_hline(yintercept = 0, linetype = 3, color="grey80")
	c <- c + theme(legend.position="right")
#  c <- c + scale_color_manual(values = c("darkgrey","black"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())
	c <- c  + ylab("difference in RA compared to healthy samples") +
		xlab("")
	#print(c)
#	dev.off()
	return(c)
}


plotDifferenceGenus <- function(df, file="figures/difference_class.pdf", jColors){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	# remove NA
	df <- df[which(!is.na(df.diff.m$g)),]
	# add type to colors
	jColors2 <- c("black", "grey")
	names(jColors2) <- c("AIH", "control")
	jColors <- c(jColors, jColors2)
	c <- ggplot(df, aes(x=reorder(g, -value), y=value*100))
 c <- c + geom_line(aes(group=g, color=c), size=2)
	c <- c + geom_point(size=5, shape=73,aes(color=variable)) + coord_flip()
	c <- c + scale_color_manual(values = jColors)
#  c <- c + scale_color_manual(values = c(control="grey", AIH="black"))
	c <- c + geom_hline(yintercept = 0, linetype = 3, color="grey80")
	c <- c + theme(legend.position="right")
#  c <- c + scale_color_manual(values = c("darkgrey","black"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())
	c <- c  + ylab("difference in RA compared to healthy samples") +
		xlab("")
	print(c)
	dev.off()
	return(c)
}

plotDifferenceSpecies <- function(df, file="figures/difference_species.pdf", jColors){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	# remove NA
	df <- df[which(!is.na(df.diff.m$g)),]
	# add type to colors
	jColors2 <- c("black", "grey")
	names(jColors2) <- c("AIH", "control")
	jColors <- c(jColors, jColors2)
	c <- ggplot(df, aes(x=reorder(s, -value), y=value*100))
 c <- c + geom_line(aes(group=s, color=c), size=2)
	c <- c + geom_point(size=5, shape=73,aes(color=variable)) + coord_flip()
	c <- c + scale_color_manual(values = jColors)
#  c <- c + scale_color_manual(values = c(control="grey", AIH="black"))
	c <- c + geom_hline(yintercept = 0, linetype = 3, color="grey80")
	c <- c + theme(legend.position="right")
#  c <- c + scale_color_manual(values = c("darkgrey","black"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())
	c <- c  + ylab("difference in RA compared to healthy samples") +
		xlab("")
	print(c)
	dev.off()
	return(c)
}

plotDifferenceOrder <- function(df, file="figures/difference_order.pdf"){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	c <- ggplot(df, aes(x=o, y=value*100,color=factor(variable)))
	c <- c + theme(legend.position="right")
	c <- c + geom_hline(yintercept = 0, linetype = 2, color="darkgreen")
	c <- c + geom_point() + coord_flip()
	c <- c + geom_line(aes(group=o), color="grey60") + scale_color_manual(values = c("darkblue","darkred"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())

	print(c)
	dev.off()
	return(c)
}


plotDifferenceFamily <- function(df, file="figures/difference_family.pdf"){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	c <- ggplot(df, aes(x=f, y=value*100,color=factor(variable)))
	c <- c + theme(legend.position="right")
	c <- c + geom_hline(yintercept = 0, linetype = 2, color="darkgreen")
	c <- c + geom_point() + coord_flip()
	c <- c + geom_line(aes(group=f), color="grey60") + scale_color_manual(values = c("darkblue","darkred"))
	c <- c + theme_bw() + theme(panel.grimajor = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())
	c <- c + scale_y_continuous(limits = c(-10, 20), breaks = c(-10,-5,0,5,10,15,20))

	print(c)
	dev.off()
	return(c)
}



plotAbundanceClass <- function(df, file="figures/abundance_class.pdf"){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	c <- ggplot(df, aes(x=c, y=value*100, fill=factor(variable)))
	c <- c + geom_bar(stat="identity", width = 0.4) + coord_flip()
#  c <- c + scale_fill_manual(values = c("darkblue","darkred"))
	c <- c + scale_fill_manual(values = c("darkgreen","darkblue", "darkred"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
																panel.grid.minor = element_blank(),
																panel.border = element_blank(),
																panel.background = element_blank())

	c <- c + theme(legend.position="right")
	print(c)
	dev.off()
	return(c)
}

plotAbundanceOrder <- function(df, file="figures/abundance_order.pdf"){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	c <- ggplot(df, aes(x=o, y=value*100, fill=factor(variable)))
	c <- c + geom_bar(stat="identity", width = 0.4) + coord_flip()
	c <- c + scale_fill_manual(values = c("darkgreen","darkblue", "darkred"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())

	c <- c + theme(legend.position="right")
	print(c)
	dev.off()
	return(c)
}

plotAbundanceFamily <- function(df, file="figures/abundance_family.pdf"){
	require(ggplot2)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	c <- ggplot(df, aes(x=f, y=value*100, fill=factor(variable)))
	c <- c + geom_bar(stat="identity", width = 0.4) + coord_flip()
	c <- c + scale_fill_manual(values = c("darkgreen","darkblue", "darkred"))
	c <- c + theme_bw() + theme(panel.grid.major = element_blank(),
															panel.grid.minor = element_blank(),
															panel.border = element_blank(),
															panel.background = element_blank())
	c <- c + theme(legend.position="right")
	print(c)
	dev.off()
	return(c)
}

# generate triplot
triplotTaxaClass <- function(this.df, file="figures/triplot_class.pdf",
														 label=TRUE, jColors){
	require(ggplot2)
	require(ggtern)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	this.df$average <- rowMeans(cbind(this.df$healthy,this.df$control,this.df$AIH))
	c <- ggtern(this.df,aes(x=healthy, y=control, z= AIH, color=c))
	c <- c + geom_point(aes(size = average)) # add point AIH
	c <- c + scale_color_manual(values = jColors)
#  c <- c + scale_color_manual("", values = this.df$col)  #+  guides(colour=FALSE)
	c <- c + guides(size = FALSE) + theme_minimal() # theme stuff   guides(colour=FALSE)
	c <- c + theme(legend.key = element_blank(),
								strip.background = element_rect(colour="white", fill="#FFFFFF") )
	if(label){
		c <- c + geom_text(aes(label = c)) # add text lable
	}
	print(c)
	dev.off()
	return(c)
}

# generate triplot
triplotTaxaPhylum <- function(this.df, file="figures/triplot_phylum.pdf",
														 label=TRUE, jColors){
	require(ggplot2)
	require(ggtern)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	this.df$average <- rowMeans(cbind(this.df$healthy,this.df$control,this.df$AIH))
	c <- ggtern(this.df,aes(x=healthy, y=control, z= AIH, color=p))
	c <- c + geom_point(aes(size = average)) # add point AIH
	c <- c + scale_color_manual(values = jColors)
#  c <- c + scale_color_manual("", values = this.df$col)  #+  guides(colour=FALSE)
	c <- c + guides(size = FALSE) + theme_minimal() # theme stuff   guides(colour=FALSE)
	c <- c + theme(legend.key = element_blank(),strip.background = element_rect(colour="white", fill="#FFFFFF") )
	if(label){
		c <- c + geom_text(aes(label = p)) # add text lable
	}
	print(c)
	dev.off()
	return(c)
}

# generate triplot
triplotTaxaGenus <- function(this.df, file="figures/triplot_class.pdf",
														 label=TRUE, jColors){
	require(ggplot2)
	require(ggtern)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	this.df$average <- rowMeans(cbind(this.df$healthy,this.df$control,this.df$AIH))
	c <- ggtern(this.df,aes(x=healthy, y=control, z= AIH, color=c))
	c <- c + geom_point(aes(size = average)) # add point AIH
	c <- c + scale_color_manual(values = jColors)
#  c <- c + scale_color_manual("", values = this.df$col)  #+  guides(colour=FALSE)
	c <- c + guides(size = FALSE) + theme_minimal() # theme stuff   guides(colour=FALSE)
	c <- c + theme(legend.key = element_blank(),
								strip.background = element_rect(colour="white", fill="#FFFFFF") )
	if(label){
		c <- c + geom_text(aes(label = g)) # add text lable
	}
	print(c)
	dev.off()
	return(c)
}

# generate triplot
triplotTaxaOrder <- function(df, file="figures/triplot_order.pdf",
														 label=TRUE,
														 alpha=1){
	require(ggplot2)
	require(ggtern)
	print(paste("plot written to", file))
	pdf(file, width = 7, height = 4)
	df$average <- rowMeans(cbind(df$healthy,df$control,df$AIH))
	c <- ggtern(df,aes(healthy,control,AIH, color=paste(p,c)))
	c <- c + geom_point(aes(size = average), alpha=alpha) # add point AIH
	c <- c + scale_color_manual("", values = df$col)
	c <- c + guides(size=FALSE) + theme_minimal() # theme stuff  + guides(colour=FALSE)
	c <- c + theme(legend.key = element_blank(),
								 strip.background = element_rect(colour="white", fill="#FFFFFF") )
	if(label){
		c <- c + geom_text(aes(label = o)) # add text lable
	}
	print(c)
	dev.off()
	return(c)
}

### Set up a blank theme
theme_none <- theme(
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	panel.background = element_blank(),
	axis.title.x = element_text(colour=NA),
	axis.title.y = element_blank(),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
	axis.line = element_blank()
	#axis.ticks.length = element_blank()
)

#' generates plot with logFC information and average abundance for species
#' and class level.
speciesFC <- function(df, order.by="none"){
	require(ggplot2)
	if ("name" %in% colnames(df)){ # there a 'names' column - use the annot.phylum column
		if (order.by == "none"){
			q <- ggplot(df, aes(x=logFC, y=name))
		} 
		if (order.by == "pvalue"){
			q <- ggplot(df, aes(x=reorder(logFC, -P.Value), y=name))
		}
		q <- q + geom_line(aes(group=name), size=2, color="grey90")
	} else { # there is no name column, use the annot.phylum column 
		if (order.by == "none"){
			q <- ggplot(df, aes(x=logFC, y=annot.phylum))	
		} 
		if (order.by == "pvalue"){
						q <- ggplot(df, aes(x=reorder(logFC, -P.Value), y=annot.phylum))
		}
		q <- q + geom_line(aes(group=annot.phylum), size=2, color="grey90")
	}
	q <- q + geom_vline(xintercept=0, color="black", size=1)
	q <- q + geom_point(aes(size=AveExpr, color=type), alpha=1)
	q <- q + scale_colour_manual(values=c("#7f7f7fff", "black"))
	q <- q + geom_text(aes(label=sig), color="red", size=7)
	q <- q + theme_minimal()		
	return(q)
}

