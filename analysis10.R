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
otu.biom <- normalizeBiomTable(otu.biom)

familyData.raw  <- aggTax(otu.biom, lvl = "X3", out = "MRexperiment", log=FALSE, norm=T)
mat.familyData.raw  <- aggTax(otu.biom, lvl = "X3", out = "matrix", log=FALSE, norm=T)

mat.raw <- as.data.frame(mat.familyData.raw)

data <- as.data.frame(scale(mat.raw, center=FALSE, scale=colSums(mat.raw)))


mapping <- read.table("data/mapping_cat_clean.txt", sep=";")

AIH_ids <- mapping[which(mapping$V3 == "AIH") ,]$V1
healthy_ids <- mapping[which(mapping$V3 == "control") ,]$V1
control_ids <- mapping[which(mapping$V3 == "others") ,]$V1

# AIH
idx <- match(AIH_ids, names(data))
aih_mean <- rowMeans(data[,idx]) 
aih_sd <- apply(data[,idx], 1, sd )

# control
idx <- match(control_ids, names(data))
control_mean <- rowMeans(data[,idx]) 
control_sd <- apply(data[,idx], 1, sd )

# healthy
idx <- match(healthy_ids, names(data))
healthy_mean <- rowMeans(data[,idx]) 
healthy_sd <- apply(data[,idx], 1, sd )

df <- data.frame(rownames(data),
                 aih = aih_mean,
                 control = control_mean,
                 healthy = healthy_mean)


df_sd <- data.frame(rownames(data),
                    aih = aih_sd,
                    control = control_sd,
                    healthy = healthy_sd)
# plot
df_mean <- df
df2 <- df
df2$rownames.data. <- NULL
df <- df[which(rowSums(df2) > 0.01),]
df_sd <- df_sd[which(rowSums(df2) > 0.01),]

df_sd <- df_sd[which(rownames(df) != "f__"),]
df <- df[which(rownames(df) != "f__"),]

df.m <- melt(df)

df.m$sd <- melt(df_sd)$value

df.m$name <- rownames(df.m)
p <- ggplot(df.m, aes(reorder(rownames.data., +value), value, fill=variable))
p <- p + geom_errorbar(aes(ymin=value, ymax=value+sd), width=.25,
                       position=position_dodge(.9))
p <- p + geom_bar(color="black", stat = "identity", position = position_dodge()) + coord_flip() + theme_classic()
p


