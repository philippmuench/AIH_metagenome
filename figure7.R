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
dist_mat <- read.table("data/OTU_table/bray_curtis_otu_table.txt")

# read mapping
mapping <- read.table("data/mapping_cat.txt", sep=";", header=F)


# get sample ID of AIH/FAP
idstoremove <- mapping[which(mapping$V2 == "HIV" |
                                    mapping$V2 == "FAP" ),]

# remove HIV, FAP
dist_mat2 <- dist_mat[-match(idstoremove$V1, colnames(dist_mat))
                      ,-match(idstoremove$V1, colnames(dist_mat))]



# get sample types from mapping file
mapping.names.ordered <- mapping[match(rownames(dist_mat2), mapping$V1),]$V1
mapping.types.ordered <- mapping[match(rownames(dist_mat2), mapping$V1),]$V3

obs.pcoa <- pcoa(as.dist(dist_mat2))

types.ordered <- mapping[match(rownames(obs.pcoa$vectors), mapping$V1),]$V3
colorvalues <- c("green", "red", "blue")
df <- data.frame(MDS1=obs.pcoa$vectors[,1], MDS2=obs.pcoa$vectors[,2], genotype=types.ordered)

#c <- drawPcoa(obs.pcoa, df, sample.names= mapping.types.ordered, colorvalues=colorvalues, label=F)

p <- ggplot(df, aes(x=MDS1,y=MDS2, label=mapping.names.ordered))
p <- p + geom_point(aes(size=4, alpha=0.8, colour= mapping.types.ordered)) # draw points

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


pdf("results/figure7/figure_7.pdf", width=6, height=4)
p
dev.off()

# add label
p <- p + geom_label(size = 2, aes(fill = mapping.types.ordered), colour = "white", fontface = "bold")

pdf("results/figure7/figure_7_label.pdf", width=6, height=4)
p
dev.off()


