# cleanup
rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggtern)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation

# path to files
path.L <- "data/otu_table_L3.txt" # path to L3 table

# load OTU table
raw.data <- read.table(path.L , header=T)

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
bar <- barplotClass(df, "results/figure2/figure_2_bar.pdf", tooSmall = 0.2)
dif <- plotDifferenceClass(df.diff.m,"results/figure2/figure_2_diff.pdf", jColors)
tri <- triplotTaxaClass(df, "results/figure2/figure_2_triplot.pdf", label=F, jColors)

# arrange all 3 in one plot
pdf("figures/class_level.pdf", width=12, height=10)
print(grid.arrange(bar,  arrangeGrob(tri, dif),ncol=2))
dev.off()
