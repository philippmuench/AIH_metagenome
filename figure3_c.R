# cleanup
rm(list=ls())

require(grid)
require(reshape)
require(reshape2)
require(gridExtra)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/heatmap.functions.R") # heatmap specific functions

# path to files
path.L <- "data/OTU_table/OTU-table-L7.txt" # path to OTU table
path.sample.matrix <- "data/mapping_cat.txt"
marker.ifapb <- read.table("data/marker/corr_iFAPB.csv", sep=";") # file with correlation information of ifapb marker to samples
marker.lbp <- read.table("data/marker/corr_LBP.csv", sep=";") # file with correlation information of lbp marker to samples
marker.scd14 <- read.table("data/marker/corr_sCD14.csv", sep=";") # file with correlation information of sCD14.csv marker to samples

# path to output files
out.heatmap <- "results/figure3/figure_3_c.pdf"

# get sample types
sample.num <- read.table(path.sample.matrix, sep=";")

# load OTU table
raw.data <- read.table(path.L , header=T)

# remove HIV and FAP
raw.data <- removeSampleFromDist(raw.data,"HIV", as.matrix(sample.num[match(colnames(raw.data),as.character(as.matrix(sample.num$V1))),]$V2))
raw.data <- removeSampleFromDist(raw.data,"FAP", as.matrix(sample.num[match(colnames(raw.data),as.character(as.matrix(sample.num$V1))),]$V2))

# basic manipulation of otu table
raw.data <- makeOtu(raw.data, method = "proportion")

# remove genera with rel ab. < 1%
raw.data <- removeRelAb(raw.data, 0.01)

# transform data for d
raw.data.t <- t(raw.data)

# process tax ID
raw.data <- processTaxaID(raw.data.t)

# remove all without species lvl annotation
raw.data.s <- raw.data[which(!is.na(raw.data$s)),]
rownames(raw.data.s) <- paste(raw.data.s$g, raw.data.s$s)

# extract taxa information for grouping
class <- raw.data.s$c
phylum <- raw.data.s$p
raw.data.s$k <- NULL; raw.data.s$p <- NULL; raw.data.s$c <- NULL;
raw.data.s$o <- NULL; raw.data.s$f <- NULL; raw.data.s$g <- NULL;
raw.data.s$s <- NULL;

# get sample type
type <- data.frame(type = sample.num[match(colnames(raw.data.s), sample.num$V1),]$V3)
type2 <- data.frame(type = sample.num[match(colnames(raw.data.s), sample.num$V1),]$V2)
mat <- as.data.frame( t(raw.data.s))

# sort by dendogram
dd.col <- as.dendrogram(hclust(dist(mat)))
col.ord <- order.dendrogram(dd.col)
dd.row <- as.dendrogram(hclust(dist(t(mat))))
row.ord <- order.dendrogram(dd.row)
xx <- as.matrix(mat[col.ord, row.ord]) # no scaling
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]

# prepare data frame
df$sample <- xx_names[[1]]
df$sample <- with(df, factor(sample, levels=sample, ordered=TRUE))
mdf <- melt(df, id.vars=c("sample"))
mdf$type <- sample.num[match(mdf$sample, sample.num$V1),]$V3

# generate ranges for better plotting colors
mdf$value1 <- cut(mdf$value,breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1),right = FALSE)

# get adiv values
#df.adiv.observed <- getAdiv(index_type="observed_otus")
#df.adiv.chao1 <- getAdiv(index_type="chao1")
#df.adiv.shannon <- getAdiv(index_type="shannon")

# order by dendogram
#df.adiv.observed <- df.adiv.observed[match(rownames(df), df.adiv.observed$sample),]
#df.adiv.chao1 <- df.adiv.chao1[match(rownames(df), df.adiv.chao1$sample),]
#df.adiv.shannon <- df.adiv.shannon[match(rownames(df), df.adiv.shannon$sample),]
#df.adiv.observed$sample <- with(df.adiv.observed, factor(sample, levels=sample, ordered=TRUE))
#df.adiv.chao1$sample <- with(df.adiv.chao1, factor(sample, levels=sample, ordered=TRUE))
#df.adiv.shannon$sample <- with(df.adiv.shannon, factor(sample, levels=sample, ordered=TRUE))

### healthy
marker.lbp$marker <- "lbp"
marker.scd14$marker <- "scd14"
marker.ifapb$marker <- "ifapb"

# get corr to marker
marker.f.aih <- processMarker(marker.ifapb, marker.scd14, marker.lbp, type="AIH")
marker.f.control <- processMarker(marker.ifapb, marker.scd14, marker.lbp, type="control")
marker.f.healthy <- processMarker(marker.ifapb, marker.scd14, marker.lbp, type="healthy")
marker.f.aih$type <- "AIH"
marker.f.control$type <- "control"
marker.f.healthy$type <- "healtyh"
marker <- rbind(marker.f.aih, marker.f.control, marker.f.healthy)
marker.m <- melt(marker)

# process adiv
#df.adiv.observed$method <- "observed"
#df.adiv.chao1$method <- "chao1"
#df.adiv.shannon$method <- "shannon"
#df.adiv.all <- rbind(df.adiv.observed, df.adiv.chao1, df.adiv.shannon)

# generate heatmap components
a <- drawHeatmapCat(mdf) # get the heatmap
c <- drawDendogram(dd.row, ddata_x, sample.num, annot=F)
d <- drawDendogramCat(dd.col, ddata_y, sample.num, annot=F)
h <- drawMarkerHeatmap(marker.m)

# save heatmap components
ggsave(a, file="results/figure3/figure_3_c_heatmap.pdf", width=6, height=8)
#ggsave(h, file="figures/heatmap/marker.pdf", width=6, height=8)
ggsave(d, file="results/figure3/figure_3_c_dendogram_col.pdf", width=6, height=8)
ggsave(c, file="results/figure3/figure_3_c_dendogram_row.pdf", width=6, height=8)
