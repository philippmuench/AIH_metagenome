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
otu.biom <- removeFeaturesbasedOnMRcounts(otu.biom, 5)

# normalize otu table based on quantile estimate
#otu.biom <- normalizeBiomTable(otu.biom)

data <- as.data.frame(returnAppropriateObj(otu.biom, TRUE, FALSE))

data <- as.data.frame(scale(data, center=FALSE, scale=colSums(data)))


bool_data <- data > 0.01
abundance_threshold <- apply(bool_data, 1, any)

data <- data[abundance_threshold,]

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
fit <- lmFit(df_aih_healthy_control, design=as.matrix(aih_design))
fit2 <- eBayes(fit)
fit.top.aih <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top.aih$star <- add.significance.stars(fit.top.aih$adj.P.Val)

# healthy 
fit <- lmFit(df_aih_healthy_control, design=as.matrix(healthy_design))
fit2 <- eBayes(fit)
fit.top.healthy <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top.healthy$star <- add.significance.stars(fit.top.healthy$adj.P.Val)

# control
fit <- lmFit(df_aih_healthy_control, design=as.matrix(control_design))
fit2 <- eBayes(fit)
fit.top.control <- topTable(fit2, number=Inf, adjust.method="fdr")
fit.top.control$star <- add.significance.stars(fit.top.control$adj.P.Val)




require(ggplot2)
require(ggtern)
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
write.table(aih_sig, file="results/figure2/aih_db.csv", sep='\t')
pie(table(as.character(as.matrix(aih_sig$class))))

# control
control_sig$logFC <- fit.top.aih[match(control_sig$rownames.data., rownames(fit.top.aih)),]$logFC
control_sig$AveExpr <- fit.top.aih[match(control_sig$rownames.data., rownames(fit.top.aih)),]$AveExpr
control_sig$t <- fit.top.aih[match(control_sig$rownames.data., rownames(fit.top.aih)),]$t
control_sig$P.Value <- fit.top.aih[match(control_sig$rownames.data., rownames(fit.top.aih)),]$P.Value
control_sig$adj.P.Val <- fit.top.aih[match(control_sig$rownames.data., rownames(fit.top.aih)),]$adj.P.Val

control_sig$species <- annot[match(control_sig$rownames.data., annot$OTU),]$X7
control_sig$genus <- annot[match(control_sig$rownames.data., annot$OTU),]$X6
control_sig$family <- annot[match(control_sig$rownames.data., annot$OTU),]$X5
control_sig$order <- annot[match(control_sig$rownames.data., annot$OTU),]$X4
control_sig$class <- annot[match(control_sig$rownames.data., annot$OTU),]$X3
control_sig$phylum <- annot[match(control_sig$rownames.data., annot$OTU),]$X2
pie(table(as.character(as.matrix(control_sig$class))))
write.table(control_sig, file="results/figure2/control_db.csv", sep='\t')

# healthy
healthy_sig$logFC <- fit.top.aih[match(healthy_sig$rownames.data., rownames(fit.top.aih)),]$logFC
healthy_sig$AveExpr <- fit.top.aih[match(healthy_sig$rownames.data., rownames(fit.top.aih)),]$AveExpr
healthy_sig$t <- fit.top.aih[match(healthy_sig$rownames.data., rownames(fit.top.aih)),]$t
healthy_sig$P.Value <- fit.top.aih[match(healthy_sig$rownames.data., rownames(fit.top.aih)),]$P.Value
healthy_sig$adj.P.Val <- fit.top.aih[match(healthy_sig$rownames.data., rownames(fit.top.aih)),]$adj.P.Val

healthy_sig$species <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X7
healthy_sig$genus <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X6
healthy_sig$family <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X5
healthy_sig$order <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X4
healthy_sig$class <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X3
healthy_sig$phylum <- annot[match(healthy_sig$rownames.data., annot$OTU),]$X2
pie(table(as.character(as.matrix(healthy_sig$class))))
write.table(healthy_sig, file="results/figure2/healthy_db.csv", sep='\t')



c <- ggtern(df,aes(x=healthy, y=control, z= AIH, color=color))
c <- c + geom_point(aes(size = average), alpha=1) # add point AIH

c <- c + guides(size = FALSE) + theme_classic() # theme stuff   guides(colour=FALSE)
c <- c + theme(legend.key = element_blank(),
               strip.background = element_rect(colour="white", fill="#FFFFFF") )

