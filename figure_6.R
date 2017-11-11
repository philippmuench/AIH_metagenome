# cleanup
rm(list=ls())

require(reshape2)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/alpha.diversity.functions.R")
source("src/correlation.functions.R")

# path to files
path.L <- "data/OTU_table/OTU-table-L7.txt" # path to L3 table

# path to file with mapping information
path.sample.matrix <- "data/mapping_cat.txt"

# path to file with mapping information
path.marker.matrix <- "data/marker/marker.csv"

# get sample types
sample.num <- read.table(path.sample.matrix, sep=";")

# get marker information
sample.marker<- read.table(path.marker.matrix, sep=";", header=T)

# load OTU table
raw.data <- read.table(path.L , header=T)

# get adiv values
shannon <- getAdiv(index_type="shannon")
chao1 <- getAdiv(index_type="chao1")
observed <- getAdiv(index_type="observed_otus")

# match adiv with makrer
shannon <- matchMarker(shannon, sample.marker)
chao1 <- matchMarker(chao1, sample.marker)
observed <- matchMarker(observed, sample.marker)

# melt for plotting
shannon.m <- melt(shannon, id.vars =c("sample", "type", "type2", "adiv"))
chao1.m <- melt(chao1, id.vars =c("sample", "type", "type2", "adiv"))
observed.m <- melt(observed, id.vars =c("sample", "type", "type2", "adiv"))

plot.shannon <- drawCorMarkerOtu(shannon.m, title="shannon", all=F)
plot.chao1 <- drawCorMarkerOtu(chao1.m, title="chao1", all=F)
plot.observed <- drawCorMarkerOtu(observed.m, title="observed", all=F)

plot.shannon.2 <- drawCorMarkerOtu(shannon.m, title="shannon", all=T)
plot.chao1.2 <- drawCorMarkerOtu(chao1.m, title="chao1", all=T)
plot.observed.2 <- drawCorMarkerOtu(observed.m, title="observed", all=T)

# print correlation values to screen
printCorrelationCoef(shannon)
printCorrelationCoef(chao1)
printCorrelationCoef(observed)

pdf(file="results/figure6/figure_6.pdf", width=4, height=8)
print(plot.shannon)
print(plot.chao1)
print(plot.observed)
dev.off()

# add "all" as lm to plot
pdf(file="results/figure6/figure_6_all.pdf", width=4, height=8)
print(plot.shannon.2)
print(plot.chao1.2)
print(plot.observed.2)
dev.off()