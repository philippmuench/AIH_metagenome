# Liver Evaluation

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/norm_libsizes_0_1.png)

## Taxonomic Bar Plots

OTU data can be further summarized and compared their abundance at different taxonomic levels based on the annotation. Users can choose to visualize the results in either area plot or bar plot (actual or proportion). OTUs without taxa designation will be collapsed into one NA (Not_Assigned) group, which could be an arbitrary mixture of OTUs across different levels.

### by Sample:

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/7_1.png)

### by Group

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/g2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/g3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/g4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/g5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/g6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/TaxBars/g7b_1.png)


## Alpha Diversity

Alpha diversity is used to measure the diversity within a sample or an ecosystem. The two most commonly used alpha-diversity measurements are Richness (count) and Evenness (distribution). MicrobiomeAnalyst provides many metrics to calculate diversity within samples. Most commonly used ones are listed below:

Observed: It estimates the amount of unique OTUs found in each sample (richness).

ACE and Chao1: They considers observed OTUs with account for unobserved species based on low-abundance OTUs (richness).

Shannon, Simpson and Fisher: These metrics accounts for both richness and evenness.

### Statistics:

Observed	p-value: 0.00020488; [T-test] statistic: 4.2005

Observed	p-value: 0.00028179; [Mann-Whitney] statistic: 576

Chao1	p-value: 0.00017764; [T-test] statistic: 4.2408

Chao1	p-value: 0.00017831; [Mann-Whitney] statistic: 575

ACE	p-value: 0.00022523; [T-test] statistic: 4.1661

ACE	p-value: 0.00030306; [Mann-Whitney] statistic: 568

Shannon	p-value: 0.00035228; [T-test] statistic: 3.8727

Shannon	p-value: 0.00035114; [Mann-Whitney] statistic: 566

Simpson	p-value: 0.00013364; [T-test] statistic: 4.1192

Simpson	p-value: 0.0004061; [Mann-Whitney] statistic: 564

Fisher	p-value: 0.00026627; [T-test] statistic: 4.1154

Fisher	p-value: 0.00022447; [Mann-Whitney] statistic: 572

### Observed

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/Observed_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/gObserved_1.png)

### Chao1

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/Chao1_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/gChao1_1.png)

### ACE

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/ACE_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/gACE_1.png)

### Shannon

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/Shannon_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/gShannon_1.png)

### Simpson

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/Simpson_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/gSimpson_1.png)

### Fisher

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/Fisher_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Alphadiversity/gFisher_1.png)


## Beta Diversity

Beta diversity represents the explicit comparison of microbial communities (in-between) based on their composition. It provides a measure of the distance or dissimilarity between each sample pair. Beta diversity is calculated for every pair of samples to generate a distance or dissimilarity matrix, reflecting the dissimilarity between certain samples. Beta diversity can be performed using the ordination based methods such as Principal Coordinates Analysis (PCoA), Nonmetric Multidimensional Scaling (NMDS) or Principal Component Analysis (PCA). Note, for CLR transformed data, PCA should be used instead of PCoA.

[PERMANOVA] R-squared: 0.046306; p-value < 0.001

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Betadiversity/PCoA_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Betadiversity/MNDS_1.png)

## Univariate Analysis

see folder:

 ![UnivariateStatisticalComparisons](/MicrobiomAnalysis/Liver-healthy-hepatopathy/UnivariateStatisticalComparisons/)

## metagenomeSeq

The algorithm is designed to determine features that are differentially abundant between two or more groups of multiple samples. It addresses the effects of both normalization and under-sampling of microbial communities on disease association detection and the testing of feature correlations. In particular, it can perform zero-inflated Gaussian fit or fitFeatureModel (for two groups only) on data after normalizing the data through cumulative sum scaling (CSS) method. For more details, please refer to the original publication: http://dx.doi.org/doi:10.1038/nmeth.2658

see folder:

 ![metagenomeSeq](/MicrobiomAnalysis/Liver-healthy-hepatopathy/metagenomeSeq/)



## DDseq2

DESeq2 are two powerful statistical methods originally developed for RNAseq data analysis. When applied to metagenomic datasets after proper data filtering and normalization, these two methods have been shown to perform equally well or better than a lot of methods designed for metagenomic datasets.

see folder

 ![DESeq2](/MicrobiomAnalysis/Liver-healthy-hepatopathy/DESeq2/)


### Examples:


## LEfSe

The LEfSe algorithm was developed for biomarker discovery and interpretation for metagenomics data. It employs Kruskal- Wallis rank sum test to detect features with significant differential abundance with regard to class labels, followed by Linear Discriminant Analysis to evaluate the relevance or effect size of differential abundant features.

See results in folder:

 ![LEfSe-Results](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/)

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/LEfSe/7_1.png)


## Random Forests

Random Forests is a powerful machine learning algorithms for classification as well as for identification of predictive features (biomakers). It operates by constructing a multitude of decision trees ("forests") at training time and predicting the class as the majority vote of the individual trees.

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c2_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c3_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c4_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c5_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c6_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c7_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f7_1.png)

### OTU

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/c8_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/RandomForests/f8_1.png)

## Core Microbiome

Core microbiome refers to the set of taxa that are detected in a high fraction of the population above a given abundance threshold. The count data is transformed to compositional (relative) abundances in order to perform such analysis.

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/CoreMicrobiome/genus_1.png)

## Hierarchical Clustering & Heatmap Visualization

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Heatmap/genus_2.png)
 
## Dendrogram (OTUs)

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/Dendrogram/OTU_1.png)

## Correlation Analysis (Genus)

 ![alt text](/MicrobiomAnalysis/Liver-healthy-hepatopathy/CorrelationAnalysis/genus_2.png)

