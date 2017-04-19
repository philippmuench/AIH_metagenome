#' by Philipp C. Münch (philipp.muench@helmholtz-hzi.de)

drawAlpha <- function(adiv, title, cat=FALSE, box=FALSE, rel=FALSE){
	mean_values <- data.frame(type = c("AIH", "control", "others"),
											value = c(mean(adiv[which(adiv$type2=="AIH"),]$value),
																mean(adiv[which(adiv$type2=="control"),]$value),
																mean(adiv[which(adiv$type2=="others"),]$value)))
	adiv$value2 <- adiv$value - mean(adiv$value)

if (rel){
	adiv$value <- adiv$value2
}

if (box){
	p <- ggplot(adiv, aes(type, value))
	p <- p + geom_boxplot() +  ggtitle(title)
	p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_minimal()

} else {
		p <- ggplot(adiv, aes(x=sample, y=value))
	p <- p + geom_bar(stat='identity') +  ggtitle(title)
	p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + mytheme

}
	if(cat){
	 p <- p + facet_grid(. ~ type, scales="free", shrink=TRUE, space="free")
	}
	#p <- p + facet_grid(type2 ~ ., scales="free", space="free", shrink =TRUE)
	#p <- p + theme_none
 # p <- p + geom_hline(data = mean_values, aes(yintercept = value, color=type)) + theme_minimal()
	#p <- p #+ coord_flip()
	return(p)
}

# boxplot by cat for heatmap figure
drawAlphaBox <- function(adiv){
	adiv <- adiv[order.dendrogram(dd.col),]
	# make sample an ordered factor
	#adiv$sample <- factor(adiv$sample, levels = adiv$sample)
 # mean_values <- data.frame(type = c("AIH", "control", "others"),
	 #                   value = c(mean(adiv[which(adiv$type2=="AIH"),]$value),
		#                            mean(adiv[which(adiv$type2=="control"),]$value),
		 #                           mean(adiv[which(adiv$type2=="others"),]$value)))
 # a <- ggplot(adiv, aes(1, value))
#  a <- a + stat_summary(geom = 'errorbar', fun.data = 'seFunc', width = 0, show_guide = F, position="dodge") + facet_grid(method ~ type, scales="free", shrink=TRUE,)
	#a <- a + stat_summary(geom = 'point', fun.y = 'mean', size = 2, shape = 21,position="dodge") + mytheme
	# a <- a + geom_hline(data = mean_values, aes(yintercept = value, color=type)) #+ theme_minimal()
	p <- ggplot(adiv, aes(1,value))
	p <- p + geom_boxplot() + facet_grid(method ~ type, scales="free", shrink=TRUE)
	#p <- p + geom_jitter()
	#p <- p + facet_grid(type2 ~ ., scales="free", space="free", shrink =TRUE)
	p <- p + theme_minimal()
	#  p <- p + geom_hline(data = mean_values, aes(yintercept = value, color=type)) + theme_minimal()
	#p <- p #+ coord_flip()
	 return(p)
}

getAdiv <- function(path.to.alpha.dif = "data/adiv.txt", path.sample.matrix = "data/mapping_cat.txt", index_type="shannon"){
	require(ggplot2)
	require(scales)
	# load alpha diversity
	adiv <- read.table(path.to.alpha.dif, header=T)
	# get sample types
	sample.num <- read.table(path.sample.matrix, sep=";")
	adiv$sample <- sample.num[match(rownames(adiv), sample.num$V1),]$V1
	adiv$type <- sample.num[match(rownames(adiv), sample.num$V1),]$V2
	adiv$type2 <- sample.num[match(rownames(adiv), sample.num$V1),]$V3
	if (index_type == "shannon"){
		df.adiv <- data.frame(sample=adiv$sample, value=adiv$shannon , type=adiv$type2, type2=adiv$type)
	}
	if (index_type == "chao1"){
		df.adiv <- data.frame(sample=adiv$sample, value=adiv$chao1 , type=adiv$type2, type2=adiv$type)
	}
	if (index_type == "observed_otus"){
		df.adiv <- data.frame(sample=adiv$sample, value=adiv$observed_otus , type=adiv$type2, type2=adiv$type)
	}
	# remove FAP and HIV samples
	df.adiv <- df.adiv[-which(df.adiv$type2=="FAP"),]
	df.adiv <- df.adiv[-which(df.adiv$type2=="HIV"),]
	df <- as.data.frame(melt(df.adiv))
	#limits <- aes(ymax = resp + se, ymin=resp - se)
	p <- ggplot(df, aes(x=sample, y=value))
	p <- p + geom_point() + facet_grid(type ~ ., scales="free", space="free", shrink =TRUE)
	p <- p + theme_minimal() + coord_flip()
	#print(p)
	# make sample an ordered factor
	df.adiv$sample <- factor(df.adiv$sample, levels = df.adiv$sample)
	#df.adiv <- df.adiv[order(match(rownames(df), df.adiv$sample)),]
	# make turkey hsd test
	aov.result <- aov(value ~ type, data=df.adiv)
 print( TukeyHSD(aov.result))
	return(df.adiv)
}

drawHeatmap <- function(df){
	require(ggplot2)
	require(RColorBrewer)
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	p <- ggplot(df, aes(x=variable, y=sample))
	p <- p + geom_tile(aes(fill=value)) + coord_flip()
	p <- p + scale_fill_gradientn(colours = myPalette(4)) + coord_flip()
	p <- p + theme(axis.line.x=element_blank(),
				axis.ticks.x=element_blank(),
				axis.text.x=element_blank(),
				axis.title.x=element_blank(),
				panel.background=element_rect(fill="white"),
				panel.grid=element_blank())
	return(p)
}

drawNothing <- function(){
	require(ggplot2)
	p <- grid.rect(gp=gpar(col="white"))
	return(p)
}

drawHeatmapCat <- function(df){
	require(ggplot2)
	require(RColorBrewer)
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	p <- ggplot(df, aes(x=sample, y=variable))
	p <- p + geom_tile(aes(fill=value)) #+ coord_flip()
	#p <- p +  scale_fill_manual(colors=myPalette(5))
	p <- p + scale_fill_gradientn(colours = myPalette(4))
	p <- p + facet_grid(.~type, space = "free", scales = "free", shrink = TRUE, drop = TRUE)
	p <- p + theme(axis.line.x=element_blank(),
	 #     axis.ticks.x=element_blank(),
		#    axis.text.x=element_blank(),
				axis.title.x=element_blank(),
				axis.text.x = element_text(angle = 90, hjust = 1),
				panel.background=element_rect(fill="white"),
				panel.grid=element_blank())
	return(p)
}

drawNothing <- function(){
	require(ggplot2)
	p <- grid.rect(gp=gpar(col="white"))
	return(p)
}

drawDendogram <- function(dd, data, sample.num, minimal=T, annot=F){
	require(ggplot2)
	require(ggdendro)
	data <- dendro_data(dd)
	p <- ggplot(segment(data))
	p <- p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
	if (minimal){
		p <- p + theme(axis.title.x=element_blank()) + theme_minimal()
		p <- p + theme(axis.line.x=element_blank(),
				axis.line.y=element_blank(),
				axis.ticks.x=element_blank(),
				axis.text.x=element_blank(),
				axis.title.x=element_blank(),
				axis.ticks.y=element_blank(),
				axis.text.y=element_blank(),
				axis.title.y=element_blank(),
				panel.background=element_rect(fill="white"),
				panel.grid=element_blank()) + coord_flip()
	}
	if (annot){
		dendo_type <- sample.num[match(data$labels$label, sample.num$V1),]$V2
		dendo_type2 <- sample.num[match(data$labels$label, sample.num$V1),]$V3
		txt_annot <- data.frame(x=data$labels$x, y=data$labels$y, type= dendo_type, type2 = dendo_type2)
	#  p <- p + geom_text(data=txt_annot, aes(x=x, y=y, label=type, hjust=0,color=type2, size=0.3, angle = 90))
		 p <- p + geom_point(data=txt_annot, aes(x=x, y=y, color=type2, size=4))
		 p <- p + geom_point(data=txt_annot, aes(x=x, y=0.05+y, color=type, size=4))
		p <- p + guides(size=FALSE)
	}
	return(p)
}

drawDendogramCat <- function(dd, data, sample.num, minimal=T, annot=F){
	require(ggplot2)
	require(ggdendro)
	data <- dendro_data(dd)
	p <- ggplot(segment(data))
	p <- p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
	if (minimal){
		p <- p + theme(axis.title.x=element_blank()) + theme_minimal()
		p <- p + theme(axis.line.x=element_blank(),
			axis.line.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.y=element_blank(),
			axis.text.y=element_blank(),
			axis.title.y=element_blank(),
			panel.background=element_rect(fill="white"),
			panel.grid=element_blank()) + coord_flip()
	}
	if (annot){
		dendo_type <- sample.num[match(data$labels$label, sample.num$V1),]$V2
		dendo_type2 <- sample.num[match(data$labels$label, sample.num$V1),]$V3
		txt_annot <- data.frame(x=data$labels$x, y=data$labels$y, type= dendo_type, type2 = dendo_type2)
		p <- p + geom_text(data=txt_annot, aes(x=x, y=y, label=type, hjust=2,color=type2, size=0.3, angle = 90))
		p <- p + geom_text(aes(data=data, x=x, y=y, label=labels$label))
		p <- p + geom_point(data=txt_annot, aes(x=x, y=y, color=type2, size=4))
		p <- p + geom_point(data=txt_annot, aes(x=x, y=0.05+y, color=type, size=4))
		p <- p + facet_grid(.~type, space = "free", scales = "free", shrink = TRUE, drop = TRUE)
		p <- p + guides(size=FALSE)
	}
	return(p)
}

#' draws corrleation of three markers to three datasets (AIH, healthy, control)
drawMarkerHeatmap <- function(df){
	require(ggplot2)
	require(RColorBrewer)
	df$value1 <- cut(df$value,breaks = c(-1,-0.7,-0.5,-0.2,0.2,0.5,0.7,1),right = FALSE)
	p <- ggplot(df, aes(x=variable, y=name))
	p <- p + geom_tile(aes(fill=value1)) + facet_grid(.~type)
	p <- p +  scale_fill_brewer(palette = "RdYlGn")
	#p <- p + scale_fill_gradientn(colours = myPalette(100))#(8))
	p <- p + theme(axis.line.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		panel.background=element_rect(fill="white"),
		panel.grid=element_blank())
	return(p)
}

processMarker <- function(marker.ifapb, marker.scd14, marker.lbp, type="AIH"){
	marker.df <- rbind(marker.ifapb, marker.scd14, marker.lbp)
	marker.df <- marker.df[which(marker.df$type==type),]
	marker <- data.frame(species = unique(marker.df$species),genus = rep("test", length(unique(marker.df$species))), cor.ifapb=rep(0, length(unique(marker.df$species))), cor.lbp=rep(0, length(unique(marker.df$species))),cor.scd14=rep(0, length(unique(marker.df$species))))
	# iterate over unique species and and use the median if there is >1 OTU for this species
	i <- 1
	marker$genus <- as.matrix(as.character(marker$genus))
	for (species in marker$species){
		subset <- marker.df[which(marker.df$species == species),]
		marker$genus[i] <- as.character(as.matrix(subset$genus[1]))
		marker[i,]$cor.ifapb <- median(subset[which(subset$marker=="ifapb"),]$cor, na.rm=T)
		marker[i,]$cor.lbp <- median(subset[which(subset$marker=="lbp"),]$cor, na.rm=T)
		marker[i,]$cor.scd14 <- median(subset[which(subset$marker=="scd14"),]$cor, na.rm=T)
		i <- i +1
	}
	# remove characters on genus and species names
	 split <- strsplit(as.character(as.matrix(marker$species)),
											";", fixed = TRUE)
	 split2 <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
	marker$species.split <- split2
	split <- strsplit(as.character(as.matrix(marker$genus)),
											";", fixed = TRUE)
	 split2 <- sapply(strsplit(sapply(split, "[", 1),"_"),"[", 3)
	marker$genus.split <- split2
	marker$name <- paste(marker$genus.split, marker$species.split)
	marker$species.split <- NULL; marker$genus.split <- NULL
	marker$species <- NULL; marker$genus <- NULL
	rightorder <- match(colnames(df), marker$name)
	marker <- marker[rightorder,]
	marker$name <- with(marker, factor(name, levels=name, ordered=TRUE))
	return(marker)
}
