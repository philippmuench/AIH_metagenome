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

familyData.raw  <- aggTax(otu.biom, lvl = "X3", out = "MRexperiment", log=FALSE)
mat.familyData.raw  <- aggTax(otu.biom, lvl = "X3", out = "matrix", log=FALSE)
# get the values for plotting
mat.raw <- as.data.frame(mat.familyData.raw)

mat.raw <- as.data.frame(scale(mat.raw, center=FALSE, scale=colSums(mat.raw)))

mat.t.raw <- as.data.frame(t(mat.raw))
#as.data.frame(t(returnAppropriateObj(biom, norm=FALSE, log=TRUE)))
#colnames(mat) <- fData(biom)$V2
mat.t.raw$type <- pData(familyData.raw)$Treatment
mat <- mat.t.raw



#mat.raw <- as.data.frame(scale(mat.raw, center=FALSE, scale=colSums(mat.raw)))
#mat.t <- as.data.frame(t(mat))
#as.data.frame(t(returnAppropriateObj(biom, norm=FALSE, log=TRUE)))
#colnames(mat) <- fData(biom)$V2
#mat.t$type <- pData(familyData)$Treatment


# statistics aih vs. healthy
mat.aih.healthy <- mat[which(pData(familyData.raw)$Treatment=="AIH" | pData(familyData.raw)$Treatment == "control"),]
design.aih.healthy <- (mat.aih.healthy$type=="AIH")*1
mat.aih.healthy <- mat.aih.healthy[,-ncol(mat.aih.healthy)]
fit <- lmFit(t(mat.aih.healthy), design=as.matrix(design.aih.healthy))
fit2 <- eBayes(fit)
fit.top <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top$star <- add.significance.stars(fit.top$adj.P.Val)
write.table(fit.top, file="results/figure11/zig_aih_vs_helathy.tsv", sep="\t", quote=F)


# statistics aih vs. control
mat.aih.control <- mat[which(pData(familyData.raw)$Treatment=="AIH" | pData(familyData.raw)$Treatment != "control"),]
design.aih.control <- (mat.aih.control$type=="AIH")*1
mat.aih.control <- mat.aih.control[,-ncol(mat.aih.control)]
fit <- lmFit(t(mat.aih.control), design=as.matrix(design.aih.control))
fit2 <- eBayes(fit)
fit.top <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top$star <- add.significance.stars(fit.top$adj.P.Val)
write.table(fit.top, file="results/figure11/zig_aih_vs_control.tsv", sep="\t", quote=F)

# statistics healthy vs. control
mat.healthy.control <- mat[which(pData(familyData.raw)$Treatment=="control" | pData(familyData.raw)$Treatment != "AIH"),]
design.healthy.control <- (mat.healthy.control$type=="control")*1
mat.healthy.control <- mat.healthy.control[,-ncol(mat.healthy.control)]
fit <- lmFit(t(mat.healthy.control), design=as.matrix(design.healthy.control))
fit2 <- eBayes(fit)
fit.top <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top$star <- add.significance.stars(fit.top$adj.P.Val)
write.table(fit.top, file="results/figure11/zig_healthy_vs_control.tsv", sep="\t", quote=F)

mat.raw[which(mat.raw$type=="control"),]$type <- "healthy"
mat.raw[which(mat.raw$type!="healthy" & mat.raw$type !="AIH"),]$type <- "control"

#mat$type <- pData(familyData)$Treatment
df <- melt(mat.raw)

#require(ggplot2)
#mytheme2 <- theme_minimal()
#mytheme2$axis.line.x <- mytheme2$axis.line.y <- mytheme2$axis.line
#a <- ggplot(df, aes(reorder(variable,-value), value, color=type))
#a <- a + stat_summary(geom = 'errorbar', fun.data = 'seFunc', width = 0, show_guide = F, position="dodge")
#a <- a + stat_summary(geom = 'point', fun.y = 'mean', size = 2, shape = 21,position="dodge") + mytheme2
#a <- a + theme(axis.text.x=element_text(angle = -45, hjust = 0))
#a <- a + scale_color_manual(values = c("black","grey50", "#00c094ff"),labels = c("AIH","control", "healthy"))+ scale_y_log10()

# plot
#pdf("figure3/genus_unpooled_log10.pdf", width=40, height=6)
#print(a)
#dev.off()
