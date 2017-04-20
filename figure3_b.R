# cleanup
rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggtern)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation

#### Generate Tri plot of genus level ####

# path to files
path.L <- "data/OTU_table/Treatment_otu_table_L6.txt" # path to L6 table

# load OTU table
raw.data <- read.table(path.L , header=T, sep=",")

# remove outgroup
raw.data$outgroup <- NULL

# basic manipulation of otu table
raw.data <- makeOtu(raw.data, method = "proportion")

# remove genera with rel ab. < 1%
raw.data <- removeRelAb(raw.data, 0.01)

# transform data for barchart
raw.data.t <- t(raw.data)

# process tax ID
df <- processTaxaID(raw.data.t)

# assign colors to taxa
df$col <-  taxCol(df$p, df$c)
jColors <- df$col
names(jColors) <- df$c
jColors <- jColors[order(names(jColors))]

# melt dataframe for plotting
df.m <- melt(df)

# calculate difference to healthy
df.diff <- processTaxaDiff(df, "healthy", c("AIH","control"))
df.diff$healthy <- NULL # remove control
df.diff.m <- melt(df.diff) # melt dataframe for plotting

# generate plots
#bar <- barplotGenus(df, "figures/genus_level_bar.pdf", tooSmall = 0.2)
#dif <- plotDifferenceGenus(df.diff.m,"results/figure3/f", jColors)
tri <- triplotTaxaGenus(df, "results/figure3/figure_3_b.pdf", label=T, jColors)

# arrange all 3 in one plot
#pdf("figures/genus_level.pdf", width=12, height=10)
#print(grid.arrange(arrangeGrob(tri, dif),ncol=1))
#dev.off()