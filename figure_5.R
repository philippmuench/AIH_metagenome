# cleanup
rm(list=ls())

require(reshape2)
require(ggplot2)


source("src/plot.functions.R") # all function for plotting
source("src/calc.functions.R") # functions for data manipulation
source("src/alpha.diversity.functions.R") # functions for alpha diversity analysis

# get values and statistics
df.observed <- getAdiv(index_type="observed_otus")
df.chao1 <- getAdiv(index_type="chao1")
df.shannon <- getAdiv(index_type="shannon")

alpha <- data.frame(sample= df.observed$sample, observed = df.observed$value,
                    chao1 = df.chao1$value, shannon = df.shannon$value,
                    cat = df.observed$type, type = df.observed$type2)

# load marker
marker <- read.table("data/parenchy.csv", header=T, sep=";")

# append mats
alpha$marker <- marker[match(alpha$sample, marker$sample), ]$paren
alpha$immuno <- marker[match(alpha$sample, marker$sample), ]$immuno

# plot
alpha.m <- melt(alpha)

a <- ggplot(alpha.m, aes(reorder(marker, - value), value, fill = cat))
a <- a + stat_summary(geom = 'errorbar', fun.data = 'seFunc', width = 0, aes(color = cat), show_guide = F)
a <- a + stat_summary(geom = 'point', fun.y = 'mean', size = 3, shape = 21) + facet_wrap(~ variable, scale="free")
a <- a + scale_fill_manual(values = c("black","#00c094ff", "grey50"),labels = c("AIH","healthy", "control"))
a <- a + scale_color_manual(values = c("black","#00c094ff", "grey50"),labels = c("AIH","healthy", "control"))
a <- a + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("results/figure5/figure_5.pdf", a)


b <- ggplot(alpha.m, aes(reorder(marker, - value), value))
b <- b + stat_summary(geom = 'errorbar', fun.data = 'seFunc', width = 0, show_guide = F)
b <- b + stat_summary(geom = 'point', fun.y = 'mean', size = 3, shape = 21) + facet_wrap(~ variable, scale="free")
b <- b + scale_fill_manual(values = c("black","#00c094ff", "grey50"),labels = c("AIH","healthy", "control"))
b <- b + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("results/figure5/figure_5_no_cat.pdf", b)

# perfrom anova and tukeyHSD
amod <- aov(chao1 ~ marker, data=alpha)
summary(amod)
coefficients(amod)
testout <- TukeyHSD(amod)
print(testout)
