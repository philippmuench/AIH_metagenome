# AIH versus healthy

 ![alt text](/MicrobiomAnalysis/AIH-healthy/norm_libsizes_0_1.png)

## Taxonomic Bar Plots

OTU data can be further summarized and compared their abundance at different taxonomic levels based on the annotation. Users can choose to visualize the results in either area plot or bar plot (actual or proportion). OTUs without taxa designation will be collapsed into one NA (Not_Assigned) group, which could be an arbitrary mixture of OTUs across different levels.

### by Sample:

### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/7_1.png)

### by Group

### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/g2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/g3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/g4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/g5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/g6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/AIH-healthy/TaxBars/g7b_1.png)


## Alpha Diversity

Alpha diversity is used to measure the diversity within a sample or an ecosystem. The two most commonly used alpha-diversity measurements are Richness (count) and Evenness (distribution). MicrobiomeAnalyst provides many metrics to calculate diversity within samples. Most commonly used ones are listed below:

Observed: It estimates the amount of unique OTUs found in each sample (richness).

ACE and Chao1: They considers observed OTUs with account for unobserved species based on low-abundance OTUs (richness).

Shannon, Simpson and Fisher: These metrics accounts for both richness and evenness.

### Statistics:

Observed	p-value: 0.028403; [T-test] statistic: -2.2923

Observed	p-value: 0.041618; [Mann-Whitney] statistic: 90

Chao1	p-value: 0.030117; [T-test] statistic: -2.2665

Chao1	p-value: 0.047883; [Mann-Whitney] statistic: 92

ACE	p-value: 0.035863; [T-test] statistic: -2.1879

ACE	p-value: 0.056186; [Mann-Whitney] statistic: 94

Shannon	p-value: 0.052649; [T-test] statistic: -2.0208

Shannon	p-value: 0.076275; [Mann-Whitney] statistic: 98

Simpson	p-value: 0.047993; [T-test] statistic: -2.1074

Simpson	p-value: 0.044123; [Mann-Whitney] statistic: 91

Fisher	p-value: 0.037503; [T-test] statistic: -2.1677

Fisher	p-value: 0.060754; [Mann-Whitney] statistic: 95


### Observed

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/Observed_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/gObserved_1.png)

### Chao1

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/Chao1_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/gChao1_1.png)

### ACE

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/ACE_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/gACE_1.png)

### Shannon

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/Shannon_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/gShannon_1.png)

### Simpson

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/Simpson_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/gSimpson_1.png)

### Fisher

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/Fisher_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Alphadiversity/gFisher_1.png)


## Beta Diversity

Beta diversity represents the explicit comparison of microbial communities (in-between) based on their composition. It provides a measure of the distance or dissimilarity between each sample pair. Beta diversity is calculated for every pair of samples to generate a distance or dissimilarity matrix, reflecting the dissimilarity between certain samples. Beta diversity can be performed using the ordination based methods such as Principal Coordinates Analysis (PCoA), Nonmetric Multidimensional Scaling (NMDS) or Principal Component Analysis (PCA). Note, for CLR transformed data, PCA should be used instead of PCoA.

[PERMANOVA] R-squared: 0.053642; p-value < 0.008

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Betadiversity/PCoA_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/Betadiversity/MNDS_1.png)

## Univariate Analysis

see folder:

 ![UnivariateStatisticalComparisons](/MicrobiomAnalysis/AIH-healthy/UnivariateStatisticalComparisons/)

| "Family"         | "Pvalues"  | "FDR"     | "Statistics" | 
|------------------|------------|-----------|--------------| 
| "Prevotellaceae" | 0.00010518 | 0.0042073 | 263          | 
 

## metagenomeSeq

The algorithm is designed to determine features that are differentially abundant between two or more groups of multiple samples. It addresses the effects of both normalization and under-sampling of microbial communities on disease association detection and the testing of feature correlations. In particular, it can perform zero-inflated Gaussian fit or fitFeatureModel (for two groups only) on data after normalizing the data through cumulative sum scaling (CSS) method. For more details, please refer to the original publication: http://dx.doi.org/doi:10.1038/nmeth.2658

see folder:

 ![metagenomeSeq](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/)

| "Order"              | "Pvalues" | "FDR"     | 
|----------------------|-----------|-----------| 
| "Desulfovibrionales" | 0.0001678 | 0.0040273 | 

| "Family"              | "Pvalues"  | "FDR"     | 
|-----------------------|------------|-----------| 
| "Desulfovibrionaceae" | 5.1672e-05 | 0.0014974 | 
| "Prevotellaceae"      | 7.3042e-05 | 0.0014974 | 

| "Genus"           | "Pvalues"  | "FDR"   | 
|-------------------|------------|---------| 
| "Aggregatibacter" | 0.00083817 | 0.03522 | 
| "Gemmiger"        | 0.0011579  | 0.03522 | 
| "Ruminococcus"    | 0.0016239  | 0.03522 | 
| "Prevotella"      | 0.0021211  | 0.03522 | 
| "Paraprevotella"  | 0.0024458  | 0.03522 | 

| "Species"                                                | "Pvalues"  | "FDR"    | 
|----------------------------------------------------------|------------|----------| 
| "Paraprevotella_clara"                                   | 0.00017653 | 0.022567 | 
| "Alistipes_indistinctus"                                 | 0.0002033  | 0.022567 | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6"             | 0.00040427 | 0.029916 | 
| "Bacteroides_plebeius"                                   | 0.00057768 | 0.032061 | 
| "OTU392-NN=Aggregatibacter_aphrophilus_EF605278-D=95_1"  | 0.00082535 | 0.036646 | 
| "Parabacteroides_goldsteinii"                            | 0.0013696  | 0.038491 | 
| "OTU810-NN=Soleaferrea_massiliensis_JX101688-D=87_6"     | 0.0013822  | 0.038491 | 
| "OTU790-NN=Barnesiella_intestinihominis_AB370251-D=84_6" | 0.0013871  | 0.038491 | 
| "Christensenella_minuta"                                 | 0.0015687  | 0.038695 | 
| "OTU1413-NN=Mesorhizobium_loti_RLU50164-D=82_1"          | 0.0020361  | 0.042041 | 
| "OTU1322-NN=Catabacter_hongkongensis_AB671763-D=85_4"    | 0.0020831  | 0.042041 | 

| "OTU"                                                    | "Pvalues"  | "FDR"      | 
|----------------------------------------------------------|------------|------------| 
| "Prevotella"                                             | 1.3856e-06 | 0.00050575 | 
| "Paraprevotella_clara"                                   | 6.3353e-06 | 0.0011562  | 
| "Streptococcus_1"                                        | 2.4898e-05 | 0.0025918  | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6_1"           | 2.8403e-05 | 0.0025918  | 
| "Alistipes_indistinctus"                                 | 0.00017398 | 0.012701   | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6"             | 0.00034135 | 0.020766   | 
| "Pasteurellaceae"                                        | 0.00040667 | 0.021205   | 
| "OTU392-NN=Aggregatibacter_aphrophilus_EF605278-D=95_1"  | 0.00053583 | 0.024447   | 
| "Christensenella_minuta"                                 | 0.0011929  | 0.041874   | 
| "OTU790-NN=Barnesiella_intestinihominis_AB370251-D=84_6" | 0.0012137  | 0.041874   | 
| "OTU1322-NN=Catabacter_hongkongensis_AB671763-D=85_4"    | 0.0014106  | 0.041874   | 
| "Parabacteroides_goldsteinii"                            | 0.0014342  | 0.041874   | 
| "Clostridiales_66"                                       | 0.0014914  | 0.041874   | 
| "Bacteroides_plebeius"                                   | 0.0016434  | 0.042846   | 
| "OTU810-NN=Soleaferrea_massiliensis_JX101688-D=87_6"     | 0.0020115  | 0.047708   | 
| "Bacteroides_1"                                          | 0.0020913  | 0.047708   | 


### Examples:

 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Aggregatibacter_1.png)

 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Deltaproteobacteria_1.png)
 
 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Desulfovibrionaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Desulfovibrionales_1.png)
 
 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Gemmiger_1.png)
 
 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Paraprevotella_1.png)

 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Prevotella_1.png)
 
 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Prevotellaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/AIH-healthy/metagenomeSeq/Ruminococcus_1.png)



## DDseq2

DESeq2 are two powerful statistical methods originally developed for RNAseq data analysis. When applied to metagenomic datasets after proper data filtering and normalization, these two methods have been shown to perform equally well or better than a lot of methods designed for metagenomic datasets.

see folder

 ![DESeq2](/MicrobiomAnalysis/AIH-healthy/DESeq2/)

| "Order"              | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|----------------------|----------|---------|------------|------------| 
| "Pasteurellales"     | -3.2828  | 0.76832 | 1.9315e-05 | 0.00046357 | 
| "Desulfovibrionales" | 2.259    | 0.74326 | 0.0023712  | 0.028454   | 

| "Family"              | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|-----------------------|----------|---------|------------|------------| 
| "Prevotellaceae"      | -4.324   | 0.67306 | 1.3241e-10 | 5.4288e-09 | 
| "Desulfovibrionaceae" | 2.6789   | 0.71741 | 0.00018838 | 0.0038617  | 
| "Pasteurellaceae"     | -2.3121  | 0.77072 | 0.0027006  | 0.032134   | 
| "unclassified"        | 1.4173   | 0.47977 | 0.0031351  | 0.032134   | 
| "Streptococcaceae"    | -1.4422  | 0.51901 | 0.0054585  | 0.04476    | 

| "Genus"          | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|------------------|----------|---------|------------|------------| 
| "Prevotella"     | -6.607   | 0.75442 | 1.9922e-18 | 1.4344e-16 | 
| "Paraprevotella" | -6.0128  | 1.0306  | 5.4077e-09 | 1.8068e-07 | 
| "Veillonella"    | -4.4314  | 0.76686 | 7.5284e-09 | 1.8068e-07 | 
| "Gemmiger"       | 2.9862   | 0.77514 | 0.00011691 | 0.0021045  | 
| "Haemophilus"    | -3.2677  | 0.86186 | 0.00014978 | 0.0021568  | 
| "Streptococcus"  | -1.6967  | 0.57683 | 0.0032681  | 0.039217   | 

| "Species"                                                 | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|-----------------------------------------------------------|----------|---------|------------|------------| 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6"              | -8.087   | 0.84046 | 6.4494e-22 | 1.4318e-19 | 
| "Paraprevotella_clara"                                    | -5.7221  | 0.99062 | 7.6371e-09 | 8.4772e-07 | 
| "OTU126-NN=Alistipes_marseilloanorexicus_JX101692-D=90_9" | 5.9543   | 1.0574  | 1.7883e-08 | 1.3233e-06 | 
| "OTU1028-NN=Ruminococcus_gnavus_JN713312-D=95_6"          | -5.9975  | 1.0785  | 2.6854e-08 | 1.4904e-06 | 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"            | -4.1668  | 0.77764 | 8.4031e-08 | 3.731e-06  | 
| "Haemophilus_parainfluenzae"                              | -3.8565  | 0.80785 | 1.8075e-06 | 5.9672e-05 | 
| "Bacteroides_plebeius"                                    | -5.4291  | 1.1392  | 1.8815e-06 | 5.9672e-05 | 
| "Streptococcus_infantis"                                  | -2.463   | 0.56944 | 1.5227e-05 | 0.00042255 | 
| "OTU1678-NN=Veillonella_dispar_GQ422726-D=94_9"           | -4.8505  | 1.2275  | 7.7614e-05 | 0.0018682  | 
| "OTU407-NN=Soleaferrea_massiliensis_JX101688-D=87_4"      | 3.5865   | 0.91207 | 8.4155e-05 | 0.0018682  | 
| "OTU1087-NN=Blautia_glucerasea_AB588023-D=93"             | 2.7722   | 0.76667 | 0.00029932 | 0.0060409  | 
| "Streptococcus_anginosus"                                 | -3.1039  | 0.89408 | 0.00051728 | 0.0090294  | 
| "Alistipes_obesi"                                         | 2.2929   | 0.6616  | 0.00052875 | 0.0090294  | 
| "OTU1033-NN=Ruminococcus_bromii_DQ882649-D=88"            | 2.3425   | 0.69611 | 0.00076491 | 0.012129   | 
| "OTU916-NN=Eubacterium_siraeum_EUBRRDV-D=93_2"            | 3.3592   | 1.0229  | 0.0010233  | 0.015145   | 
| "OTU790-NN=Barnesiella_intestinihominis_AB370251-D=84_6"  | -3.3292  | 1.0732  | 0.0019212  | 0.026365   | 
| "OTU1196-NN=Christensenella_minuta_AB490809-D=86_5"       | 3.4589   | 1.1203  | 0.002019   | 0.026365   | 
| "OTU260-NN=Clostridium_clostridioforme_AY169422-D=94_1"   | 2.0581   | 0.70641 | 0.0035742  | 0.044082   | 

| "OTU"                                                        | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|--------------------------------------------------------------|----------|---------|------------|------------| 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6"                 | -8.2212  | 0.86825 | 2.834e-21  | 1.0344e-18 | 
| "Streptococcus_1"                                            | -4.4991  | 0.65354 | 5.8084e-12 | 1.06e-09   | 
| "Paraprevotella_clara"                                       | -6.4187  | 0.99044 | 9.1314e-11 | 1.111e-08  | 
| "OTU126-NN=Alistipes_marseilloanorexicus_JX101692-D=90_9"    | 6.2502   | 1.0689  | 4.9973e-09 | 4.56e-07   | 
| "Prevotella"                                                 | -6.5325  | 1.1586  | 1.7162e-08 | 1.2528e-06 | 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"               | -4.22    | 0.75951 | 2.7579e-08 | 1.5258e-06 | 
| "OTU1028-NN=Ruminococcus_gnavus_JN713312-D=95_6"             | -6.0777  | 1.0959  | 2.9262e-08 | 1.5258e-06 | 
| "Streptococcus_infantis"                                     | -3.0751  | 0.59175 | 2.0285e-07 | 9.2551e-06 | 
| "Bacteroides_plebeius"                                       | -5.7624  | 1.1678  | 8.0366e-07 | 3.2593e-05 | 
| "Haemophilus_parainfluenzae"                                 | -3.7846  | 0.81116 | 3.0755e-06 | 0.00011225 | 
| "OTU293-NN=Clostridium_glycyrrhizinilyticum_AB233029-D=91_2" | 3.5985   | 0.80652 | 8.127e-06  | 0.00026967 | 
| "Streptococcus_anginosus"                                    | -3.9171  | 0.92916 | 2.4899e-05 | 0.00075734 | 
| "OTU407-NN=Soleaferrea_massiliensis_JX101688-D=87_4"         | 3.6656   | 0.93263 | 8.4811e-05 | 0.0023812  | 
| "OTU1678-NN=Veillonella_dispar_GQ422726-D=94_9"              | -4.7314  | 1.2387  | 0.00013367 | 0.0033148  | 
| "Gemmiger"                                                   | 3.0028   | 0.78711 | 0.00013623 | 0.0033148  | 
| "Streptococcus"                                              | -2.1998  | 0.58139 | 0.00015449 | 0.0035243  | 
| "Lachnospiraceae"                                            | 2.1302   | 0.56831 | 0.00017802 | 0.0038221  | 
| "OTU1087-NN=Blautia_glucerasea_AB588023-D=93"                | 2.8103   | 0.76659 | 0.00024636 | 0.0049956  | 
| "Bacteroides_1"                                              | -3.8807  | 1.0657  | 0.00027104 | 0.0052069  | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6_1"               | -4.4399  | 1.2471  | 0.00037078 | 0.0067667  | 
| "Eubacterium_1"                                              | 1.9365   | 0.55602 | 0.00049627 | 0.0086256  | 
| "OTU1033-NN=Ruminococcus_bromii_DQ882649-D=88"               | 2.3703   | 0.70178 | 0.0007312  | 0.012131   | 
| "OTU916-NN=Eubacterium_siraeum_EUBRRDV-D=93_2"               | 3.4714   | 1.0446  | 0.00089058 | 0.01358    | 
| "OTU790-NN=Barnesiella_intestinihominis_AB370251-D=84_6"     | -3.5869  | 1.0797  | 0.0008929  | 0.01358    | 
| "Alistipes_obesi"                                            | 2.2048   | 0.66938 | 0.00098833 | 0.01443    | 
| "OTU1196-NN=Christensenella_minuta_AB490809-D=86_5"          | 3.6108   | 1.1479  | 0.0016581  | 0.023277   | 
| "unclassified_16"                                            | 3.2389   | 1.0357  | 0.0017641  | 0.023847   | 
| "Paraprevotella_clara_1"                                     | -3.9043  | 1.2643  | 0.002015   | 0.026267   | 
| "OTU392-NN=Aggregatibacter_aphrophilus_EF605278-D=95_1"      | -3.5151  | 1.1496  | 0.0022299  | 0.028067   | 
| "OTU84-NN=Robinsoniella_peoriensis_AF445258-D=93_1"          | 2.4398   | 0.80249 | 0.0023631  | 0.02875    | 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89_2"      | 2.8165   | 0.95335 | 0.0031338  | 0.036898   | 
| "Ruminococcaceae_18"                                         | -3.1503  | 1.0941  | 0.003985   | 0.044549   | 
| "Fusicatenibacter_saccharivorans"                            | 1.8857   | 0.65568 | 0.0040277  | 0.044549   | 


## LEfSe

The LEfSe algorithm was developed for biomarker discovery and interpretation for metagenomics data. It employs Kruskal- Wallis rank sum test to detect features with significant differential abundance with regard to class labels, followed by Linear Discriminant Analysis to evaluate the relevance or effect size of differential abundant features.

See results in folder:

 ![LEfSe-Results](/MicrobiomAnalysis/AIH-healthy/LEfSe/)

| "Family"         | "Pvalues"  | "FDR"     | "AIH"            | "healthy"        | "LDAscore" | 
|------------------|------------|-----------|------------------|------------------|------------| 
| "Prevotellaceae" | 0.00023731 | 0.0097298 | 1296470.47603501 | 114263.998730267 | -5.77      | 


### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-healthy/LEfSe/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-healthy/LEfSe/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-healthy/LEfSe/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/AIH-healthy/LEfSe/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-healthy/LEfSe/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/AIH-healthy/LEfSe/7_1.png)


## Random Forests

Random Forests is a powerful machine learning algorithms for classification as well as for identification of predictive features (biomakers). It operates by constructing a multitude of decision trees ("forests") at training time and predicting the class as the majority vote of the individual trees.

### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c2_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c3_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c4_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c5_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c6_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c7_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f7_1.png)

### OTU

 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/c8_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-healthy/RandomForests/f8_1.png)

## Core Microbiome

Core microbiome refers to the set of taxa that are detected in a high fraction of the population above a given abundance threshold. The count data is transformed to compositional (relative) abundances in order to perform such analysis.

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-healthy/CoreMicrobiome/genus_1.png)

## Hierarchical Clustering & Heatmap Visualization

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Heatmap/genus_2.png)
 
## Dendrogram (OTUs)

 ![alt text](/MicrobiomAnalysis/AIH-healthy/Dendrogram/OTU_1.png)

## Correlation Analysis (Genus)

 ![alt text](/MicrobiomAnalysis/AIH-healthy/CorrelationAnalysis/genus_2.png)

