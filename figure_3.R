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

source("src/plot.functions.R")
source("src/calc.functions.R")
source("src/zig.functions.R")

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

# remove features that are not present in > 20 samples
otu.biom <- removeFeaturesbasedOnMRcounts(otu.biom, 20)

# normalize otu table based on quantile estimate
#otu.biom <- normalizeBiomTable(otu.biom)

data <- as.data.frame(returnAppropriateObj(otu.biom, norm=TRUE, log=TRUE, sl=1))

#data <- as.data.frame(scale(data, center=FALSE, scale=colSums(data)))

#bool_data <- data > 0.5
#abundance_threshold <- apply(bool_data, 1, any)
#data <- data[abundance_threshold,]

mapping <- read.table("data/mapping_cat_clean.txt", sep=";")

AIH_ids <- mapping[which(mapping$V3 == "AIH") ,]$V1
healthy_ids <- mapping[which(mapping$V3 == "control") ,]$V1
control_ids <- mapping[which(mapping$V3 == "others") ,]$V1


# AIH
idx_aih <- match(AIH_ids, names(data))
aih_mean <- rowMeans(data[,idx_aih]) 
aih_sd <- apply(data[,idx_aih], 1, sd )

# control
idx_control <- match(control_ids, names(data))
control_mean <- rowMeans(data[,idx_control]) 
control_sd <- apply(data[,idx_control], 1, sd )

# healthy
idx_healthy <- match(healthy_ids, names(data))
healthy_mean <- rowMeans(data[,idx_healthy]) 
healthy_sd <- apply(data[,idx_healthy], 1, sd )

aih_healthy_control <- c(idx_aih, idx_healthy, idx_control)

aih_design <- c(rep(1,length(idx_aih)), rep(0, length(idx_healthy)), rep(0, length(idx_control)))
healthy_design <- c(rep(0,length(idx_aih)), rep(1, length(idx_healthy)), rep(0, length(idx_control)))
control_design <- c(rep(0,length(idx_aih)), rep(0, length(idx_healthy)), rep(1, length(idx_control)))


#aih_control <- c(idx_aih, idx_control)
#aih_control_design <-  c(rep(0,length(idx_aih)), rep(1, length(idx_control)))

#healthy_control <- c(idx_healthy, idx_control)
#healthy_control_design <-  c(rep(0,length(idx_healthy)), rep(1, length(idx_control)))


df <- data.frame(rownames(data),
                 AIH = aih_mean,
                 control = control_mean,
                 healthy = healthy_mean)

#df_aih_healthy <- data[,aih_healthy]
#df_aih_control <- data[,aih_control]
#df_healthy_control <- data[,healthy_control]
df_aih_healthy_control <- data[, aih_healthy_control]

# aih vs healthy
#fit <- lmFit(df_aih_healthy, design=as.matrix(aih_healthy_design))
#fit2 <- eBayes(fit)
#fit.top_aih_healthy <- topTable(fit2, number=Inf, adjust.method="fdr")
#fit.top_aih_healthy$star <- add.significance.stars(fit.top_aih_healthy$adj.P.Val)

# aih vs control
#fit <- lmFit(df_aih_control, design=as.matrix(aih_control_design))
#fit2 <- eBayes(fit)
#fit.top_aih_control <- topTable(fit2, number=Inf, adjust.method="fdr")
#fit.top_aih_control$star <- add.significance.stars(fit.top_aih_control$adj.P.Val)

# healthy vs control
#fit <- lmFit(df_healthy_control, design=as.matrix(healthy_control_design))
#fit2 <- eBayes(fit)
#fit.top_healthy_control <- topTable(fit2, number=Inf, adjust.method="fdr")
#fit.top_healthy_control$star <- add.significance.stars(fit.top_healthy_control$adj.P.Val)

# aih 


#df_aih_healthy_control <- log2(df_aih_healthy_control)

groups<-as.factor(aih_design)
design_aih<-model.matrix(~groups)
fit <- lmFit(df_aih_healthy_control, design=design_aih)
fit2 <- eBayes(fit)
fit.top.aih <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top.aih$star <- add.significance.stars(fit.top.aih$adj.P.Val)
fit.top.aih <- fit.top.aih[which(fit.top.aih$logFC > 0),]
fit.top.aih <- fit.top.aih[which(fit.top.aih$P.Value < 0.05),]


# healthy 
groups<-as.factor(healthy_design)
design_healthy<-model.matrix(~groups)
fit <- lmFit(df_aih_healthy_control, design=design_healthy)
fit2 <- eBayes(fit)
fit.top.healthy <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top.healthy$star <- add.significance.stars(fit.top.healthy$adj.P.Val)
fit.top.healthy <- fit.top.healthy[which(fit.top.healthy$logFC > 0),]
fit.top.healthy <- fit.top.healthy[which(fit.top.healthy$P.Value < 0.05),]


# control
groups<-as.factor(control_design)
design_control<-model.matrix(~groups)
fit <- lmFit(df_aih_healthy_control, design=design_control)
fit2 <- eBayes(fit)
fit.top.control <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top.control$star <- add.significance.stars(fit.top.control$adj.P.Val)
fit.top.control <- fit.top.control[which(fit.top.control$logFC > 0),]
fit.top.control <- fit.top.control[which(fit.top.control$P.Value < 0.05),]


require(ggplot2)
require(ggtern)
df_zig <- df
df$average <- rowMeans(cbind(df$healthy,df$control,df$AIH))
df$color <- "normal"

aih_ones <- which(!is.na(match(rownames(df),rownames(fit.top.aih[which(fit.top.aih$P.Value < 0.05),]))))
control_ones <- which(!is.na(match(rownames(df),rownames(fit.top.control[which(fit.top.control$P.Value < 0.05),]))))
healthy_ones <- which(!is.na(match(rownames(df),rownames(fit.top.healthy[which(fit.top.healthy$P.Value < 0.05),]))))

df[aih_ones,]$color <- "aih"
df[control_ones,]$color <- "control"
df[healthy_ones,]$color <- "healthy"

aih_sig <- df[aih_ones,]
control_sig <- df[control_ones,]
healthy_sig <- df[healthy_ones,]


annot <- as.data.frame(fData(otu.biom))
annot$OTU <- rownames(annot)

#annotate fit table
## AIH
aih_sig$logFC <- fit.top.aih[match(aih_sig$rownames.data., rownames(fit.top.aih)),]$logFC
aih_sig$AveExpr <- fit.top.aih[match(aih_sig$rownames.data., rownames(fit.top.aih)),]$AveExpr
aih_sig$t <- fit.top.aih[match(aih_sig$rownames.data., rownames(fit.top.aih)),]$t
aih_sig$P.Value <- fit.top.aih[match(aih_sig$rownames.data., rownames(fit.top.aih)),]$P.Value
aih_sig$adj.P.Val <- fit.top.aih[match(aih_sig$rownames.data., rownames(fit.top.aih)),]$adj.P.Val

aih_sig$species <- annot[match(aih_sig$rownames.data., annot$OTU),]$X7
aih_sig$genus <- annot[match(aih_sig$rownames.data., annot$OTU),]$X6
aih_sig$family <- annot[match(aih_sig$rownames.data., annot$OTU),]$X5
aih_sig$order <- annot[match(aih_sig$rownames.data., annot$OTU),]$X4
aih_sig$class <- annot[match(aih_sig$rownames.data., annot$OTU),]$X3
aih_sig$phylum <- annot[match(aih_sig$rownames.data., annot$OTU),]$X2
write.table(aih_sig, file="results/figure_3_aih_db.csv", sep='\t')
pie(table(as.character(as.matrix(aih_sig$class))))

# control
control_sig$logFC <- fit.top.control[match(control_sig$rownames.data., rownames(fit.top.control)),]$logFC
control_sig$AveExpr <- fit.top.control[match(control_sig$rownames.data., rownames(fit.top.control)),]$AveExpr
control_sig$t <- fit.top.control[match(control_sig$rownames.data., rownames(fit.top.control)),]$t
control_sig$P.Value <- fit.top.control[match(control_sig$rownames.data., rownames(fit.top.control)),]$P.Value
control_sig$adj.P.Val <- fit.top.control[match(control_sig$rownames.data., rownames(fit.top.control)),]$adj.P.Val

control_sig$species <- annot[match(control_sig$rownames.data., annot$OTU),]$X7
control_sig$genus <- annot[match(control_sig$rownames.data., annot$OTU),]$X6
control_sig$family <- annot[match(control_sig$rownames.data., annot$OTU),]$X5
control_sig$order <- annot[match(control_sig$rownames.data., annot$OTU),]$X4
control_sig$class <- annot[match(control_sig$rownames.data., annot$OTU),]$X3
control_sig$phylum <- annot[match(control_sig$rownames.data., annot$OTU),]$X2
pie(table(as.character(as.matrix(control_sig$class))))
write.table(control_sig, file="results/figure_3_control_db.csv", sep='\t')

# healthy
healthy_sig$logFC <- fit.top.healthy[match(healthy_sig$rownames.data., rownames(fit.top.healthy)),]$logFC
healthy_sig$AveExpr <- fit.top.healthy[match(healthy_sig$rownames.data., rownames(fit.top.healthy)),]$AveExpr
healthy_sig$t <- fit.top.healthy[match(healthy_sig$rownames.data., rownames(fit.top.healthy)),]$t
healthy_sig$P.Value <- fit.top.healthy[match(healthy_sig$rownames.data., rownames(fit.top.healthy)),]$P.Value
healthy_sig$adj.P.Val <- fit.top.healthy[match(healthy_sig$rownames.data., rownames(fit.top.healthy)),]$adj.P.Val

healthy_sig$species <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X7
healthy_sig$genus <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X6
healthy_sig$family <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X5
healthy_sig$order <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X4
healthy_sig$class <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X3
healthy_sig$phylum <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X2
pie(table(as.character(as.matrix(healthy_sig$class))))
write.table(healthy_sig, file="results/figure_3_healthy_db.csv", sep='\t')

df_clean <- df[which(df$color != "normal"),]


c <- ggtern(df_clean,aes(x=healthy, y=control, z= AIH, color=color))
c <- c + geom_point(aes(size = average), alpha=1) # add point AIH

c <- c + guides(size = FALSE) + theme_classic() # theme stuff   guides(colour=FALSE)
c <- c + theme(legend.key = element_blank(),
               strip.background = element_rect(colour="white", fill="#FFFFFF") )
c

ggsave(filename = "results/figure_3_triplot.pdf", c)

