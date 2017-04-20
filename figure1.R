# cleanup
rm(list=ls())

require(ggplot2)
require(reshape2)

source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/alpha.diversity.functions.R") # functions for alpha diversity analysis

# get values and statistics
df.observed <- getAdiv(index_type="observed_otus")
df.chao1 <- getAdiv(index_type="chao1")
df.shannon <- getAdiv(index_type="shannon")

# generate the overview plots
plot.observed <- plotAdiv(df.observed)
plot.chao1 <- plotAdiv(df.chao1)
plot.shannon <- plotAdiv(df.shannon)

# save as pdf
pdf(file="results/figure1/figure_1_index_observed.pdf", width=8, height=6)
print(plot.observed)
dev.off()
pdf(file="results/figure1/figure1_index_chao1.pdf", width=8, height=6)
print(plot.chao1)
dev.off()
pdf(file="results/figure1/figure_1_index_shannon.pdf", width=8, height=6)
print(plot.shannon)-
dev.off()
