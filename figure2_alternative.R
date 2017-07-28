# cleanup
rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggtern)
library(matrixStats)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation

# path to files
path.L <- "data/OTU_table/final_L3.txt" # path to L3 table

mapping <- read.table("data/mapping_cat_clean.txt", sep=";")

AIH_ids <- mapping[which(mapping$V3 == "AIH") ,]$V1
healthy_ids <- mapping[which(mapping$V3 == "control") ,]$V1
control_ids <- mapping[which(mapping$V3 == "others") ,]$V1

# load OTU table
data <- as.data.frame(read.table(path.L , header=T, sep = ","))

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

df <- data.frame(data$ID,
                 aih = aih_mean,
                 control = control_mean,
                 healthy = healthy_mean)


df_sd <- data.frame(data$ID,
                 aih = aih_sd,
                 control = control_sd,
                 healthy = healthy_sd)
# plot
df_mean <- df
rownames(df_mean) <- df$data.ID
df_mean$data.ID <- NULL
df <- df[which(rowSums(df_mean) > 0.05),]
df_sd <- df_sd[which(rowSums(df_mean) > 0.05),]
df.m <- melt(df)

df.m$sd <- melt(df_sd)$value


p <- ggplot(df.m, aes(reorder(data.ID, +value), value, fill=variable))
p <- p + geom_errorbar(aes(ymin=value, ymax=value+sd), width=.25,
                       position=position_dodge(.9))
p <- p + geom_bar(color="black", stat = "identity", position = position_dodge()) + coord_flip() + theme_classic()
p
