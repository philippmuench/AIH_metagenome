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

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/zig.functions.R") # function for Zero-infated Gaussian mixture model

# path to files
path.L <- "data/OTU_table/final_L3.txt"
path.sample.matrix <- "data/mapping.txt"

# load OTU table
raw.data <- read.table(path.L , header=T, sep=",")

# basic manipulation of otu table
raw.data <- makeOtu(raw.data, method = "proportion")

# remove genera with rel ab. < 1%
#raw.data <- removeRelAb(raw.data, 0.01)

# transform data for barchart
raw.data.t <- t(raw.data)

clin <- load_phenoData(path.sample.matrix,tran = TRUE, sep = ";")
ord <- match(colnames(raw.data.t),rownames(clin))
clin <- clin[ord, ]
phenotypeData <- AnnotatedDataFrame(clin)

taxa <- read.delim("data/OTU_table/L3_taxa.txt",stringsAsFactors = FALSE, sep=",", header=F)
rownames(taxa) <- taxa$V1
OTUdata <- AnnotatedDataFrame(taxa)
biom <- newMRexperiment(raw.data.t, phenoData=phenotypeData, featureData=OTUdata, normFactors=NULL)

biom <- removeSamplefromBiom(biom, "HIV")
biom <- removeSamplefromBiom(biom, "FAP")

# remove features that are not present in > 5 samples
#biom <- removeFeaturesbasedOnMRcounts(biom, 5)

# normalize otu table based on quantile estimate
biom <- normalizeBiomTable(biom)

#### ZIG analysis ####
# this part will output two tables to the ../results/ folder

# get model
zig.treatment <- getZigModBasedOnGrouping(biom,pData(biom)$Treatment)
zig.treatment.fit <- zig.treatment$fit
zig.treatment.design <- zig.treatment$fit$design

# define what to compare
contrast.1 <- makeContrasts(groupingAIH - groupingcontrol, levels = zig.treatment.design) # AIH vs. healthy

# get all significant ones
top.contrast.1 <- fitModelWithContrast(biom, zig.treatment.fit, contrast.1, onlysig = F, phylum=T)
write.table(top.contrast.1, file="results/figure2/figure_2_c_aih_vs_healthy.csv", quote=F, sep=";")

# define what to compare
contrast.2<- makeContrasts(groupingAIH - (groupingALD + groupingalkohol + groupingHCV+ groupingNASH + groupingothers)/5, levels = zig.treatment.design)

# get all significant ones
top.contrast.2 <- fitModelWithContrast(biom, zig.treatment.fit, contrast.2, onlysig = F, phylum=T)
write.table(top.contrast.2, file="results/figure2/figure_2_c_aih_vs_control.csv", quote=F, sep=";")

# generate class index for plotting
classIndex = list(AIH = which(pData(biom)$Treatment == "AIH"))
classIndex$healthy =which(pData(biom)$Treatment == "control")
classIndex$others = which(pData(biom)$Treatment != "control" & pData(biom)$Treatment != "AIH" )

# generate volcano plot
#which(top.contrast.1$adj.P.Val < 0.001)

top.contrast.1$logP <- log10(top.contrast.1$adj.P.Val)
top.contrast.1.best <- top.contrast.1[which(top.contrast.1$logP < -3),]

vol.healthy <- createVolcan(top.contrast.1, phylum=T, avExprAsSize=T)
vol.control <- createVolcan(top.contrast.2, phylum=T, avExprAsSize=T)

#pdf("results/figure2/zig_class_volcano.pdf", width = 10, height = 5)
#grid.arrange(vol.healthy, vol.control, ncol=2)
#dev.off()