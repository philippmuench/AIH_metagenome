# function to load the alpha diversity
getAdiv <- function(path.to.alpha.dif = "data/alpha_diversity/alpha_diversity.txt", path.sample.matrix = "data/sample_mapping_cat_adiv.txt", index_type="shannon"){
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

  # make sample an ordered factor
  df.adiv$sample <- factor(df.adiv$sample, levels = df.adiv$sample)
  #df.adiv <- df.adiv[order(match(rownames(df), df.adiv$sample)),]
  # make turkey hsd test
  aov.result <- aov(value ~ type, data=df.adiv)
 print( TukeyHSD(aov.result))
  return(df.adiv)
}

# function to plot alpha diversity
plotAdiv <- function(df){
	require(ggplot2)
	a <- ggplot(df, aes(reorder(type2, -value), value, fill = type))
	a <- a + stat_summary(geom = 'errorbar', fun.data = 'seFunc', width = 0, aes(color = type), show_guide = F)
	a <- a + stat_summary(geom = 'point', fun.y = 'mean', size = 3, shape = 21)
	a <- a + scale_fill_manual(values = c("black","#00c094ff", "grey50"),labels = c("AIH","healthy", "control"))
	a <- a + scale_color_manual(values = c("black","#00c094ff", "grey50"),labels = c("AIH","healthy", "control"))
	a <- a + theme_minimal()
	return(a)
}
