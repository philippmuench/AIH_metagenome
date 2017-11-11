# cleanup
rm(list=ls())

require(ggplot2)
require(reshape2)

source("src/plot.functions.R")
source("src/calc.functions.R") 
source("src/alpha.diversity.functions.R")

# get values and statistics
df.observed <- getAdiv(index_type="observed_otus")
df.chao1 <- getAdiv(index_type="chao1")
df.shannon <- getAdiv(index_type="shannon")

# generate the overview plots
plot.observed <- plotAdiv(df.observed)
plot.chao1 <- plotAdiv(df.chao1)
plot.shannon <- plotAdiv(df.shannon)

pdf(file="results/figure_1c_chao1.pdf", width=3, height=2)
print(plot.chao1)
dev.off()

# to plot also shannon index and oberserved OTUs uncommend these lines
# save as pdf
#pdf(file="results/figure_1c_observed.pdf", width=3, height=2)
#print(plot.observed)
#dev.off()

#pdf(file="results/figure_1c_shannon.pdf", width=3, height=2)
#print(plot.shannon)
#dev.off()
