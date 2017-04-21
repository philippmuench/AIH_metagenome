# code adapted from http://enterotype.embl.de/enterotypes.html

# cleanup
rm(list=ls())

# load libs
require(phyloseq)
require(cluster)
require(clusterSim)
require(ade4)
require(ggplot2)
require(RColorBrewer)
require(ape)
require(ellipse)

source("src/pcoa.functions.R")

# load distance matrix
dist <- read.table("data/OTU_table/bray_curtis_otu_table.txt")

# read mapping
mapping <- read.table("data/mapping_cat.txt", sep=";", header=F)

# get sample types from mapping file
mapping.names.ordered <- mapping[match(rownames(dist), mapping$V1),]$V1
mapping.types.ordered <- mapping[match(rownames(dist), mapping$V1),]$V3

obs.pcoa <- pcoa(as.dist(dist))

types.ordered <- mapping[match(rownames(obs.pcoa$vectors), mapping$V1),]$V3
colorvalues <- c("green", "red", "blue")
df <- data.frame(MDS1=obs.pcoa$vectors[,1], MDS2=obs.pcoa$vectors[,2], genotype=types.ordered)


c <- drawPcoa(obs.pcoa, df, sample.names= mapping.types.ordered, colorvalues=colorvalues)


pdf("results/figure7/figure_7.pdf", width=6, height=4)
c
dev.off()

