# cleanup
rm(list=ls())

require(limma)
require(metagenomeSeq)
require(biom)
require(phyloseq)
require(ggtern)
require(reshape2)
require(colorspace)
require(grid)
require(gridExtra)
require(ggrepel)
require(pander)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/zig.functions.R") # function for Zero-infated Gaussian mixture model

# path to files
path.biom <- "data/zig_genus/otu.final.biom" # path to .biom table with taxonomic annotations in json format

#### prepare input file ####

# read in biom files
otu.biom <- load_biom(path.biom)

# normalize otu table based on quantile estimate
otu.biom <- normalizeBiomTable(otu.biom)

# remove HIV and FAP from biom file
otu.biom <- removeSamplefromBiom(otu.biom, "HIV")
otu.biom <- removeSamplefromBiom(otu.biom, "FAP")

# remove features that are not present in > 10 samples
#otu.biom <- removeFeaturesbasedOnMRcounts(otu.biom, 10)


# normalize otu table based on quantile estimate
#otu.biom <- normalizeBiomTable(otu.biom)

familyData.raw  <- aggTax(otu.biom, lvl = "X6", out = "MRexperiment", log=TRUE)
mat.familyData.raw  <- aggTax(otu.biom, lvl = "X6", out = "matrix", log=TRUE)
# get the values for plotting
mat.raw <- as.data.frame(mat.familyData.raw)
mat.raw <-as.data.frame(scale(mat.raw, center=FALSE, scale=colSums(mat.raw)))
mat.t.raw <- as.data.frame(t(mat.raw))
#as.data.frame(t(returnAppropriateObj(biom, norm=FALSE, log=TRUE)))
#colnames(mat) <- fData(biom)$V2
mat.t.raw$type <- pData(familyData.raw)$Treatment
mat.raw <- mat.t.raw

familyData  <- aggTax(otu.biom, lvl = "X6", out = "MRexperiment", log=TRUE)
mat.familyData  <- aggTax(otu.biom, lvl = "X6", out = "matrix", log=TRUE)

# get the values for plotting
mat <- as.data.frame(mat.familyData)
#mat <- mat/colSums(mat)
mat <- as.data.frame(scale(mat, center=FALSE, scale=colSums(mat)))
mat.t <- as.data.frame(t(mat))
#as.data.frame(t(returnAppropriateObj(biom, norm=FALSE, log=TRUE)))
#colnames(mat) <- fData(biom)$V2
mat.t$type <- pData(familyData)$Treatment

mat <- mat.t
threshold <- 0.5
df <- as.data.frame(mat)
df2 <- df
df2$type <- NULL
df2 <- df2[,which(colSums(df2) > threshold)]
df2$type <- df$type

df3 <- df2
df2$g__ <- NULL
df2$no_match <- NULL
dfx <- df2
dfx$type <- NULL
dfx2 <- dfx[, which(colSums(dfx) > 2)]
dfx2_b <- dfx[, which(colSums(dfx) < 2)]

dfx2$type <- df2$type
dfx2_b$type <- df2$type


df2.m <- melt(dfx2)
df2.m[which(df2.m$type == "control"),]$type <- "healthy"
df2.m[which(df2.m$type != "healthy" & df2.m$type != "AIH"),]$type <- "control"

df3.m <- melt(dfx2_b)
df3.m[which(df3.m$type == "control"),]$type <- "healthy"
df3.m[which(df3.m$type != "healthy" & df3.m$type != "AIH"),]$type <- "control"

#first
p1 <- ggplot(df2.m, aes(reorder(variable, value), value, color=type ))
p1 <- p1 + geom_boxplot(width = 0.4,  outlier.size=1) + coord_flip() + theme_classic()
#p1 <- p1 + scale_y_log10()
p1 <- p1 +  scale_color_manual(values=c("#000000", "#7f7f7f", "#00c094"))
p1 <- p1 + labs(x = "class level", y="relative abundance") 

#second
p2 <- ggplot(df3.m, aes(reorder(variable, value), value, color=type ))
p2 <- p2 + geom_boxplot(width = 0.8,  outlier.size=1) + coord_flip() + theme_classic()
#p2 <- p2 + scale_y_log10()
p2 <- p2 +  scale_color_manual(values=c("#000000", "#7f7f7f", "#00c094"))
p2 <- p2 + labs(x = "class level", y="relative abundance") 

## convert plots to gtable objects
library(gtable)
library(grid) # low-level grid functions are required
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size="first") # stack the two plots
grid.newpage()
grid.draw(g)

ggsave(filename = "results/figure_s3.pdf")

