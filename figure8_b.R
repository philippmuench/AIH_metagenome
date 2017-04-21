# cleanup
rm(list=ls())

# load packages
require(vegan)
require(pander)
require(calibrate)
require(ggplot2)

source("src/constrained.pcoa.functions.R")

# load distance matrix
dist_mat <- read.table("data/OTU_table/bray_curtis_otu_table.txt")

# load mapping file
sample.types <- read.table("data/mapping_cat.txt", sep=";", header=F)

# get sample ID of AIH/FAP
idstoremove <- sample.types[which(sample.types$V2 == "HIV" |
                     sample.types$V2 == "FAP" ),]

# remove HIV, FAP
dist_mat2 <- dist_mat[-match(idstoremove$V1, colnames(dist_mat))
         ,-match(idstoremove$V1, colnames(dist_mat))]

# load mapping file
#sample.types <- read.table("data/mapping_cat.txt", sep=";", header=F)
sample.types <- read.table("data/parenchy.csv", sep=";", header=T)

d <- data.frame(Compartment=sample.types$paren,
                Experiment=sample.types$paren,
                Description=sample.types$immuno,
                row.names=sample.types$sample)

# load grouping
sample.names <- rownames(dist_mat2)
sample.names.ordered <- sample.names[match(colnames(dist_mat2),
                                           as.character(as.matrix(sample.names)))]
sample.types.ordered <- sample.types[match(rownames(dist_mat2),
                                           sample.types$sample ),]$paren
sample.sites.ordered <- sample.types[match(colnames(dist_mat2),
                                           as.character(as.matrix(sample.types$sample))),]$paren

sample.sites.ordered.cex <- rep(0, length(sample.sites.ordered))
sample.sites.ordered.cex[which(sample.sites.ordered == "mit_fibrose")] <- 1
sample.sites.ordered.cex[which(sample.sites.ordered == "mit_zirrhose")] <- 2
sample.sites.ordered.cex[which(sample.sites.ordered == "keine_hepatopathie")] <- 3
sample.sites.ordered.cex[which(sample.sites.ordered == "hepatop_ohne_veraenderung")] <- 4


# run capscale
cap <- capscale(as.dist(dist_mat2) ~ sample.sites.ordered, add=T, sqrt.dist=T)
permanova <- anova.cca(cap)


# generate variability tables and calculate confidence intervals for the variance
cap.tbl <- variability_table(cap)
cap.var <- cap_var_props(cap)
cap.ci <-  cca_ci(cap) # get confidence interval vrom CAA object
print(cap.var)
print(cap.ci)

cap.wa <- cap$CCA$wa # extract the weighted average (sample) scores
# extract centroids of constrained factor
cap.centroids <- cap$CCA$centroids[, 1:2]


# prepare df
df <- as.data.frame(cbind(sample=rownames(dist_mat2),cap$CCA$wa[,1:2])) # plot pcoa
constraint <- "Bray-Curtis - constrained by liver status"

type.names  <- sample.sites.ordered
colorvalues <- c("green", "red", "blue","grey")

# plot
g <- ggplot(df, aes(x=as.numeric(as.matrix(CAP1)),
                    y=as.numeric(as.matrix(CAP2)), label=sample), color=type.names)
g <- g + geom_point(size= 4, alpha = 0.8, aes(color= type.names))
g <- g + scale_color_manual(values = colorvalues)
g <- g + theme_bw() + geom_hline(yintercept = 0,linetype = 3)
g <- g + geom_vline(xintercept = 0,linetype = 3)
g <- g + ggtitle(paste(constraint, ": [", format(cap.tbl["constrained", "proportion"]*100, digits=2),
                       "% of variance; \n P < ", format(permanova[1,4], digits=2),
                       "; 95% CI = ", format(cap.ci[1]*100, digits=2),
                       "%, ", format(cap.ci[2]*100, digits=2), "%]", sep="") )
g <- g + xlab(paste("Constrained PCoA 1 (", format(cap.var[1]*100, digits=4), " %)", sep=""))
g <- g + ylab(paste("Constrained PCoA 2 (", format(cap.var[2]*100, digits=4), " %)", sep=""))
g <- g + theme(plot.title = element_text(size = 8),
               axis.title.x = element_text(size = 8),
               axis.title.y = element_text(size = 8))
g


g <- g + theme(legend.position="none") + coord_equal()
pdf("results/figure8/figure_8_b_print.pdf", width=6, height=4)
g
dev.off()

g <- g + theme(legend.position="right") + coord_equal()
pdf("results/figure8/figure_8_constrained_by_liver.pdf", width=6, height=4)
g
dev.off()

# add label
g <- g + geom_label(size = 2, aes(fill = type.names), colour = "white", fontface = "bold")
pdf("results/figure8/figure_8_constrained_by_liver_label.pdf", width=6, height=4)
g
dev.off()

