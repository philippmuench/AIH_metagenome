# Liver Evaluation

 ![alt text](/MicrobiomAnalysis/Lever-all/norm_libsizes_0_1.png)

## Taxonomic Bar Plots

by Sample:

### Phylum

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/7_1.png)

by Group

### Phylum

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/g2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/g3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/g4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/g5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/g6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Lever-all/TaxBars/g7_1.png)


## Alpha Diversity

Alpha diversity is used to measure the diversity within a sample or an ecosystem. The two most commonly used alpha-diversity measurements are Richness (count) and Evenness (distribution). MicrobiomeAnalyst provides many metrics to calculate diversity within samples. Most commonly used ones are listed below:

Observed: It estimates the amount of unique OTUs found in each sample (richness).

ACE and Chao1: They considers observed OTUs with account for unobserved species based on low-abundance OTUs (richness).

Shannon, Simpson and Fisher: These metrics accounts for both richness and evenness.

Statistics:

Observed	p-value: 4.9197e-05; [ANOVA] F-value: 9.2855

Observed	p-value: 0.00040517; [Kruskal-Wallis] statistic: 18.173

Chao1	p-value: 2.1717e-05; [ANOVA] F-value: 10.158

Chao1	p-value: 0.00017699; [Kruskal-Wallis] statistic: 19.912

ACE	p-value: 2.4285e-05; [ANOVA] F-value: 10.037

ACE	p-value: 0.00022085; [Kruskal-Wallis] statistic: 19.448

Shannon	p-value: 0.0084092; [ANOVA] F-value: 4.3277

Shannon	p-value: 0.0057115; [Kruskal-Wallis] statistic: 12.552

Simpson	p-value: 0.0089386; [ANOVA] F-value: 4.2735

Simpson	p-value: 0.0069664; [Kruskal-Wallis] statistic: 12.125

Fisher	p-value: 3.2094e-05; [ANOVA] F-value: 9.738

Fisher	p-value: 0.00020323; [Kruskal-Wallis] statistic: 19.623

### Observed

 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/Observed_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/gObserved_1.png)

### Chao1

 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/Chao1_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/gChao1_1.png)

### ACE

 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/ACE_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/gACE_1.png)

### Shannon

 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/Shannon_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/gShannon_1.png)

### Simpson

 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/Simpson_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/gSimpson_1.png)

### Fisher

 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/Fisher_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Alphadiversity/gFisher_1.png)


## Beta Diversity

Beta diversity represents the explicit comparison of microbial communities (in-between) based on their composition. It provides a measure of the distance or dissimilarity between each sample pair. Beta diversity is calculated for every pair of samples to generate a distance or dissimilarity matrix, reflecting the dissimilarity between certain samples. Beta diversity can be performed using the ordination based methods such as Principal Coordinates Analysis (PCoA), Nonmetric Multidimensional Scaling (NMDS) or Principal Component Analysis (PCA). Note, for CLR transformed data, PCA should be used instead of PCoA.

[PERMANOVA] R-squared: 0.10529; p-value < 0.001

 ![alt text](/MicrobiomAnalysis/Lever-all/Betadiversity/PCoA_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/Betadiversity/MNDS_1.png)

## Univariate Analysis

 ![include](/MicrobiomAnalysis/Lever-all/UnivariateStatisticalComparisons/3.csv)
 
see folder

## metagenomeSeq

The algorithm is designed to determine features that are differentially abundant between two or more groups of multiple samples. It addresses the effects of both normalization and under-sampling of microbial communities on disease association detection and the testing of feature correlations. In particular, it can perform zero-inflated Gaussian fit or fitFeatureModel (for two groups only) on data after normalizing the data through cumulative sum scaling (CSS) method. For more details, please refer to the original publication: 

see folder:

 ![alt text](/MicrobiomAnalysis/Lever-all/metagenomeSeq/

## DDseq2

DESeq2 are two powerful statistical methods originally developed for RNAseq data analysis. When applied to metagenomic datasets after proper data filtering and normalization, these two methods have been shown to perform equally well or better than a lot of methods designed for metagenomic datasets.

see folder

 ![alt text](/MicrobiomAnalysis/Lever-all/DESeq2/

## LEfSe

The LEfSe algorithm was developed for biomarker discovery and interpretation for metagenomics data. It employs Kruskal- Wallis rank sum test to detect features with significant differential abundance with regard to class labels, followed by Linear Discriminant Analysis to evaluate the relevance or effect size of differential abundant features.

See results in folder:

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/)

### Phylum

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Lever-all/LEfSe/7_1.png)

## Random Forests

Random Forests is a powerful machine learning algorithms for classification as well as for identification of predictive features (biomakers). It operates by constructing a multitude of decision trees ("forests") at training time and predicting the class as the majority vote of the individual trees.

### Phylum

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c2_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c3_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c4_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c5_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c6_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c7_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f7_1.png)

### OTU

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c8_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f8_1.png)

## Random Forests

Random Forests is a powerful machine learning algorithms for classification as well as for identification of predictive features (biomakers). It operates by constructing a multitude of decision trees ("forests") at training time and predicting the class as the majority vote of the individual trees.

### Phylum

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c2_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c3_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c4_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c5_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c6_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c7_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f7_1.png)

### OTU

 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/c8_1.png)
 ![alt text](/MicrobiomAnalysis/Lever-all/RandomForests/f8_1.png)

## Core Microbiome

Core microbiome refers to the set of taxa that are detected in a high fraction of the population above a given abundance threshold. The count data is transformed to compositional (relative) abundances in order to perform such analysis.

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/CoreMicrobiome/genus_1.png)

## Hierarchical Clustering & Heatmap Visualization

### Genus

 ![alt text](/MicrobiomAnalysis/Lever-all/Heatmap/genus_2.png)
 
## Dendrogram (OTUs)

 ![alt text](/MicrobiomAnalysis/Lever-all/Dendrogram/OTU_1.png)

## Correlation Analysis (Genus)

 ![alt text](/MicrobiomAnalysis/Lever-all/CorrelationAnalysis/genus_2.png)
