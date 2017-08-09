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

otu.biom = filterData(otu.biom, present = 58, depth = 10)


# remove HIV and FAP from biom file
otu.biom <- removeSamplefromBiom(otu.biom, "HIV")
otu.biom <- removeSamplefromBiom(otu.biom, "FAP")

# normalize otu table based on quantile estimate
otu.biom <- normalizeBiomTable(otu.biom)

data <- as.data.frame(returnAppropriateObj(otu.biom, TRUE, FALSE))

data <- as.data.frame(scale(data, center=FALSE, scale=colSums(data)))

#bool_data <- data > 0.01
#abundance_threshold <- apply(bool_data, 1, any)

#data <- data[abundance_threshold,]

mapping <- read.table("data/mapping_cat_clean.txt", sep=";")

AIH_ids <- mapping[which(mapping$V3 == "AIH") ,]$V1
healthy_ids <- mapping[which(mapping$V3 == "control") ,]$V1
control_ids <- mapping[which(mapping$V3 == "others") ,]$V1


# AIH
idx_aih <- match(AIH_ids, names(data))
aih_mean <- rowMeans(data[,idx_aih]) 
aih_min <- apply(data[,idx_aih], 1, min )
aih_max <- apply(data[,idx_aih], 1, max )


# control
idx_control <- match(control_ids, names(data))
control_mean <- rowMeans(data[,idx_control]) 
control_min <- apply(data[,idx_control], 1, min )
control_max <- apply(data[,idx_control], 1, max )


# healthy
idx_healthy <- match(healthy_ids, names(data))
healthy_mean <- rowMeans(data[,idx_healthy]) 
healthy_min <- apply(data[,idx_healthy], 1, min )
healthy_max <- apply(data[,idx_healthy], 1, max )

aih_healthy_control <- c(idx_aih, idx_healthy, idx_control)


df <- data.frame(rownames(data),
                 AIH = aih_mean,
                 control = control_mean,
                 healthy = healthy_mean,
                 AIH_min=aih_min,
                 AIH_max=aih_max,
                 control_min = control_min,
                 control_max = control_max,
                 healthy_min = healthy_min,
                 healthy_max = healthy_max
                 )


settings = zigControl(maxit = 10, verbose = FALSE)

dat = filterData(otu.biom, present = 20, depth = 10)
pmat <- pData(dat)
pmat[which(pmat$Treatment == "control"),]$Treatment <- "healthy"
pmat[which(pmat$Treatment != "healthy" & pmat$Treatment != "AIH"),]$Treatment <- "control"
mod <- model.matrix( ~ Treatment, data=pmat) 
colnames(mod) <- c("AIH", "control", "healthy")
res = fitZig(obj = dat, mod = mod, settings = settings)
zigFit = res$fit
finalMod = res$fit$design


# healthy
contrast.matrix = makeContrasts(healthy - (control + AIH)/2, AIH - (control + healthy)/2, control - (AIH + healthy)/2,  levels = finalMod)
healthy_fit = contrasts.fit(zigFit, contrast.matrix)
healthy_fit = eBayes(healthy_fit)


zig.top.healthy <- topTable(healthy_fit, number=10000, coef=1)

#AIH
zig.top.aih <- topTable(healthy_fit, number=10000, coef=2)


#control
zig.top.control <- topTable(healthy_fit, number=10000, coef=3)



df_zig <- df

# process zig
df_zig$average <- rowMeans(cbind(df_zig$healthy,df_zig$control,df_zig$AIH))
df_zig$color <- "normal"


aih_mark <- rownames(zig.top.aih[which(zig.top.aih$P.Value < 0.01),])
control_mark <- rownames(zig.top.control[which(zig.top.control$P.Value < 0.01),])
healthy_mark <- rownames(zig.top.healthy[which(zig.top.healthy$P.Value < 0.01),])

df_zig[match(aih_mark, rownames(df_zig)),]$color <- "AIH"
df_zig[match(healthy_mark, rownames(df_zig)),]$color <- "healthy"
df_zig[match(control_mark, rownames(df_zig)),]$color <- "control"

#aih_sig_zig <- df_zig[aih_ones_zig,]
#control_sig_zig <- df_zig[control_ones_zig,]
#healthy_sig_zig <- df_zig[healthy_ones_zig,]

# see http://www.ggtern.com/d/2.2.0/geom_errorbarX.html
annot <- as.data.frame(fData(otu.biom))
annot$OTU <- rownames(annot)

df_zig_clean <- df_zig[which(df_zig$color != "normal"),]

d <- ggtern(df_zig_clean,aes(x=healthy, y=control, z= AIH, color=color))
d <- d + geom_point(aes(size = average), alpha=1) # add point AIH
d <- d + geom_errorbarR(aes(Rmin = AIH_min , Rmax = AIH_max, width =  0.01), color = 'darkred') 

d <- d + guides(size = FALSE) + theme_classic() # theme stuff   guides(colour=FALSE)
d <- d + theme(legend.key = element_blank(),
               strip.background = element_rect(colour="white", fill="#FFFFFF") )
d

