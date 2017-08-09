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

require(ggplot2)
require(ggtern)
source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/zig.functions.R") # function for Zero-infated Gaussian mixture model

# path to files
path.biom <- "data/zig_genus/otu.final.biom" # path to .biom table with taxonomic annotations in json format

#### prepare input file ####

# read in biom files
otu.biom <- load_biom(path.biom)

# remove HIV and FAP from biom file
otu.biom <- removeSamplefromBiom(otu.biom, "HIV")
otu.biom <- removeSamplefromBiom(otu.biom, "FAP")

# normalize otu table based on quantile estimate
otu.biom <- normalizeBiomTable(otu.biom)

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


df <- data.frame(rownames(data),
                 AIH = aih_mean,
                 control = control_mean,
                 healthy = healthy_mean)



# healthy
data_healthy = filterData(otu.biom, present = 10, depth = 1)
data_healthy <- cumNorm(data_healthy, p = 0.5)
pd_healthy <- pData(data_healthy)
pd_healthy[which(pd_healthy$Treatment == "control"),]$Treatment <- "group2"
pd_healthy[which(pd_healthy$Treatment != "group2"),]$Treatment <- "group1"
mod_healthy <- model.matrix(~1 + Treatment, data = pd_healthy) # todo: binary model matrix
data_healthy = fitFeatureModel(data_healthy, mod_healthy, B = 10)
zig.top.healthy <- MRcoefs(data_healthy, number = 10000)


#AIH
data_aih = filterData(otu.biom, present = 10, depth = 1)
data_aih <- cumNorm(data_aih, p = 0.5)
pd_aih <- pData(data_aih)
pd_aih[which(pd_aih$Treatment == "AIH"),]$Treatment <- "case"
pd_aih[which(pd_aih$Treatment != "case"),]$Treatment <- "control"
mod_aih <- model.matrix(~1 + Treatment, data = pd_aih) # todo: binary model matrix
data_aih = fitFeatureModel(data_aih, mod_aih,  B = 10)
zig.top.aih <- MRcoefs(data_aih, number = 10000)

#control
data_control = filterData(otu.biom, present = 10, depth = 1)
data_control <- cumNorm(data_control, p = 0.5)
pd_control <- pData(data_control)
pd_control[which(pd_control$Treatment != "control" & pd_control$Treatment != "AIH"),]$Treatment <- "case"
pd_control[which(pd_control$Treatment != "case"),]$Treatment <- "control"
mod_control <- model.matrix(~1 + Treatment, data = pd_control) # todo: binary model matrix
data_control = fitFeatureModel(data_control, mod_control)
zig.top.control <- MRcoefs(data_control, number = 10000)


df_zig <- df

# process zig
df_zig$average <- rowMeans(cbind(df_zig$healthy,df_zig$control,df_zig$AIH))
df_zig$color <- "normal"


aih_mark <- rownames(zig.top.aih[which(zig.top.aih$pvalues < 0.05),])
#control_mark <- rownames(zig.top.control[which(zig.top.control$pvalues < 0.05),])
#healthy_mark <- rownames(zig.top.healthy[which(zig.top.healthy$pvalues < 0.05),])

df_zig[match(aih_mark, rownames(df_zig)),]$color <- "AIH"
#df_zig[match(healthy_mark, rownames(df_zig)),]$color <- "healthy"
#df_zig[match(control_mark, rownames(df_zig)),]$color <- "control"

#aih_sig_zig <- df_zig[aih_ones_zig,]
#control_sig_zig <- df_zig[control_ones_zig,]
#healthy_sig_zig <- df_zig[healthy_ones_zig,]


annot <- as.data.frame(fData(otu.biom))
annot$OTU <- rownames(annot)

df_zig_clean <- df_zig[which(df_zig$color != "normal"),]

d <- ggtern(df_zig_clean,aes(x=healthy, y=control, z= AIH, color=color))
d <- d + geom_point(aes(size = average), alpha=1) # add point AIH

d <- d + guides(size = FALSE) + theme_classic() # theme stuff   guides(colour=FALSE)
d <- d + theme(legend.key = element_blank(),
               strip.background = element_rect(colour="white", fill="#FFFFFF") )
d

