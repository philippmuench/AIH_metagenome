# Liver Evaluation

 ![alt text](/MicrobiomAnalysis/Liver-all/norm_libsizes_0_1.png)

## Taxonomic Bar Plots

OTU data can be further summarized and compared their abundance at different taxonomic levels based on the annotation. Users can choose to visualize the results in either area plot or bar plot (actual or proportion). OTUs without taxa designation will be collapsed into one NA (Not_Assigned) group, which could be an arbitrary mixture of OTUs across different levels.

### by Sample:

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/7_1.png)

### by Group

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/g2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/g3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/g4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/g5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/g6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-all/TaxBars/g7b_1.png)


## Alpha Diversity

Alpha diversity is used to measure the diversity within a sample or an ecosystem. The two most commonly used alpha-diversity measurements are Richness (count) and Evenness (distribution). MicrobiomeAnalyst provides many metrics to calculate diversity within samples. Most commonly used ones are listed below:

Observed: It estimates the amount of unique OTUs found in each sample (richness).

ACE and Chao1: They considers observed OTUs with account for unobserved species based on low-abundance OTUs (richness).

Shannon, Simpson and Fisher: These metrics accounts for both richness and evenness.

### Statistics:

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

 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/Observed_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/gObserved_1.png)

### Chao1

 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/Chao1_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/gChao1_1.png)

### ACE

 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/ACE_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/gACE_1.png)

### Shannon

 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/Shannon_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/gShannon_1.png)

### Simpson

 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/Simpson_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/gSimpson_1.png)

### Fisher

 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/Fisher_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Alphadiversity/gFisher_1.png)


## Beta Diversity

Beta diversity represents the explicit comparison of microbial communities (in-between) based on their composition. It provides a measure of the distance or dissimilarity between each sample pair. Beta diversity is calculated for every pair of samples to generate a distance or dissimilarity matrix, reflecting the dissimilarity between certain samples. Beta diversity can be performed using the ordination based methods such as Principal Coordinates Analysis (PCoA), Nonmetric Multidimensional Scaling (NMDS) or Principal Component Analysis (PCA). Note, for CLR transformed data, PCA should be used instead of PCoA.

[PERMANOVA] R-squared: 0.10529; p-value < 0.001

 ![alt text](/MicrobiomAnalysis/Liver-all/Betadiversity/PCoA_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/Betadiversity/MNDS_1.png)

## Univariate Analysis

see folder:

 ![UnivariateStatisticalComparisons](/MicrobiomAnalysis/Liver-all/UnivariateStatisticalComparisons/)

| "Class"               | "Pvalues" | "FDR"   | "Statistics" | 
|-----------------------|-----------|---------|--------------| 
| "Deltaproteobacteria" | 0.0012969 | 0.01686 | 15.715       | 

| "Order"              | "Pvalues" | "FDR"    | "Statistics" | 
|----------------------|-----------|----------|--------------| 
| "Desulfovibrionales" | 0.0012969 | 0.027235 | 15.715       | 

| "Species"                                                    | "Pvalues"  | "FDR"     | "Statistics" | 
|--------------------------------------------------------------|------------|-----------|--------------| 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89"        | 2.3411e-05 | 0.0044012 | 24.135       | 
| "OTU506-NN=Christensenella_minuta_AB490809-D=83_8"           | 6.2399e-05 | 0.0051714 | 22.093       | 
| "OTU916-NN=Eubacterium_siraeum_EUBRRDV-D=93_2"               | 8.2522e-05 | 0.0051714 | 21.509       | 
| "OTU1033-NN=Ruminococcus_bromii_DQ882649-D=88"               | 0.00024542 | 0.010694  | 19.227       | 
| "OTU383-NN=Clostridium_populeti_X71853-D=91_5"               | 0.00028442 | 0.010694  | 18.917       | 
| "OTU1335-NN=Oscillibacter_valericigenes_AB238598-D=95_6"     | 0.00038077 | 0.011931  | 18.304       | 
| "Dialister_pneumosintes"                                     | 0.00060686 | 0.014621  | 17.322       | 
| "OTU661-NN=Blautia_wexlerae_EF036467-D=94_9"                 | 0.00062217 | 0.014621  | 17.269       | 
| "OTU694-NN=Clostridium_clariflavum_NR_102987_1-D=85_2"       | 0.001108   | 0.022905  | 16.049       | 
| "OTU293-NN=Clostridium_glycyrrhizinilyticum_AB233029-D=91_2" | 0.0012471  | 0.022905  | 15.799       | 
| "OTU1689-NN=Soleaferrea_massiliensis_JX101688-D=88_9"        | 0.0013402  | 0.022905  | 15.646       | 
| "OTU342-NN=Clostridium_clariflavum_NR_102987_1-D=77_3"       | 0.0016346  | 0.025609  | 15.224       | 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"               | 0.0022356  | 0.032331  | 14.558       | 
| "OTU1339-NN=Clostridium_clariflavum_NR_102987_1-D=76_8"      | 0.0024735  | 0.033216  | 14.343       | 
| "OTU1240-NN=OscillospiraClostridium_viride_X81125-D=86"      | 0.0029348  | 0.035741  | 13.978       | 
| "Parabacteroides_goldsteinii"                                | 0.0030418  | 0.035741  | 13.902       | 
| "OTU902-NN=Prevotella_loescheii_AB547688-D=88_8"             | 0.0035049  | 0.03876   | 13.599       | 

| "OTU"                                                        | "Pvalues"  | "FDR"     | "Statistics" | 
|--------------------------------------------------------------|------------|-----------|--------------| 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89_2"      | 3.6063e-05 | 0.0082247 | 23.236       | 
| "OTU506-NN=Christensenella_minuta_AB490809-D=83_8"           | 6.2399e-05 | 0.0082247 | 22.093       | 
| "OTU916-NN=Eubacterium_siraeum_EUBRRDV-D=93_2"               | 8.2522e-05 | 0.0082247 | 21.509       | 
| "unclassified_39"                                            | 0.00014119 | 0.00982   | 20.386       | 
| "Clostridiales_54"                                           | 0.00016421 | 0.00982   | 20.069       | 
| "OTU1033-NN=Ruminococcus_bromii_DQ882649-D=88"               | 0.00024542 | 0.012149  | 19.227       | 
| "OTU383-NN=Clostridium_populeti_X71853-D=91_5"               | 0.00028442 | 0.012149  | 18.917       | 
| "unclassified_4"                                             | 0.00034055 | 0.012728  | 18.538       | 
| "Bacteroidales_5"                                            | 0.00057757 | 0.016912  | 17.426       | 
| "Dialister_pneumosintes"                                     | 0.00060686 | 0.016912  | 17.322       | 
| "OTU661-NN=Blautia_wexlerae_EF036467-D=94_9"                 | 0.00062217 | 0.016912  | 17.269       | 
| "Oscillospira_4"                                             | 0.00083863 | 0.019179  | 16.638       | 
| "Lachnoclostridium_1"                                        | 0.00088375 | 0.019179  | 16.528       | 
| "Lachnospiraceae_7"                                          | 0.00089802 | 0.019179  | 16.494       | 
| "OTU694-NN=Clostridium_clariflavum_NR_102987_1-D=85_2"       | 0.001108   | 0.022087  | 16.049       | 
| "OTU293-NN=Clostridium_glycyrrhizinilyticum_AB233029-D=91_2" | 0.0012471  | 0.022262  | 15.799       | 
| "Barnesiella_10"                                             | 0.0012821  | 0.022262  | 15.74        | 
| "OTU1689-NN=Soleaferrea_massiliensis_JX101688-D=88_9"        | 0.0013402  | 0.022262  | 15.646       | 
| "Bacteria_5"                                                 | 0.0015054  | 0.022628  | 15.399       | 
| "OTU1335-NN=Oscillibacter_valericigenes_AB238598-D=95_6"     | 0.0015136  | 0.022628  | 15.388       | 
| "OTU342-NN=Clostridium_clariflavum_NR_102987_1-D=77_3"       | 0.0016346  | 0.023274  | 15.224       | 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"               | 0.0022356  | 0.030384  | 14.558       | 
| "OTU1339-NN=Clostridium_clariflavum_NR_102987_1-D=76_8"      | 0.0024735  | 0.032156  | 14.343       | 
| "OTU1240-NN=OscillospiraClostridium_viride_X81125-D=86"      | 0.0029348  | 0.034981  | 13.978       | 
| "unclassified_38"                                            | 0.00303    | 0.034981  | 13.91        | 
| "Parabacteroides_goldsteinii"                                | 0.0030418  | 0.034981  | 13.902       | 
| "OTU902-NN=Prevotella_loescheii_AB547688-D=88_8"             | 0.0035049  | 0.038814  | 13.599       | 
| "Streptococcus_1"                                            | 0.0044657  | 0.047687  | 13.081       | 

## metagenomeSeq

The algorithm is designed to determine features that are differentially abundant between two or more groups of multiple samples. It addresses the effects of both normalization and under-sampling of microbial communities on disease association detection and the testing of feature correlations. In particular, it can perform zero-inflated Gaussian fit or fitFeatureModel (for two groups only) on data after normalizing the data through cumulative sum scaling (CSS) method. For more details, please refer to the original publication: 

see folder:

 ![metagenomeSeq](/MicrobiomAnalysis/Liver-all/metagenomeSeq/)

| "Familiy"        | "Pvalues" | "FDR"    | 
|------------------|-----------|----------| 
| "Prevotellaceae" | 0.0010435 | 0.041738 | 

| "OTU"                                                   | "Pvalues"  | "FDR"     | 
|---------------------------------------------------------|------------|-----------| 
| "Prevotella"                                            | 2.2598e-05 | 0.0067568 | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6_3"          | 7.9304e-05 | 0.011856  | 
| "Bacteroidales_5"                                       | 0.00013051 | 0.013007  | 
| "Barnesiella_10"                                        | 0.00021425 | 0.016015  | 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89_2" | 0.00068136 | 0.040745  | 
| "Alistipes_indistinctus"                                | 0.00083117 | 0.04142   | 


## DDseq2

DESeq2 are two powerful statistical methods originally developed for RNAseq data analysis. When applied to metagenomic datasets after proper data filtering and normalization, these two methods have been shown to perform equally well or better than a lot of methods designed for metagenomic datasets.

see folder

 ![DESeq2](/MicrobiomAnalysis/Liver-all/DESeq2/)

| "Class"        | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"     | 
|----------------|----------|---------|------------|-----------| 
| "Fusobacteria" | 3.2854   | 0.99116 | 0.00091734 | 0.0073387 | 

| "Order"               | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|-----------------------|----------|---------|------------|------------| 
| "Gammaproteobacteria" | 2.8556   | 0.62546 | 4.9821e-06 | 3.8269e-05 | 
| "Deltaproteobacteria" | -3.2654  | 0.7183  | 5.467e-06  | 3.8269e-05 | 
| "Bacilli"             | 1.7399   | 0.44973 | 0.00010939 | 0.00041449 | 
| "Negativicutes"       | 1.9798   | 0.51431 | 0.00011843 | 0.00041449 | 
| "Fusobacteriia"       | 2.9351   | 1.1849  | 0.013245   | 0.037087   | 

| "Order"              | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|----------------------|----------|---------|------------|------------| 
| "Desulfovibrionales" | -3.4511  | 0.74176 | 3.2788e-06 | 7.2134e-05 | 
| "Pasteurellales"     | 3.6212   | 0.8155  | 8.9752e-06 | 9.8727e-05 | 
| "Lactobacillales"    | 2.0417   | 0.49013 | 3.1054e-05 | 0.00022773 | 
| "Selenomonadales"    | 1.798    | 0.53784 | 0.00082886 | 0.0045587  | 
| "Fusobacteriales"    | 3.2105   | 1.028   | 0.0017892  | 0.0078724  | 
| "Enterobacteriales"  | 2.1255   | 0.84115 | 0.011508   | 0.037611   | 
| "Bifidobacteriales"  | 1.5343   | 0.61052 | 0.011967   | 0.037611   | 

| "Family"              | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|-----------------------|----------|---------|------------|------------| 
| "Desulfovibrionaceae" | -4.3162  | 0.70563 | 9.5523e-10 | 3.8209e-08 | 
| "Streptococcaceae"    | 3.3468   | 0.59883 | 2.2849e-08 | 4.5698e-07 | 
| "Veillonellaceae"     | 3.5854   | 0.70513 | 3.6822e-07 | 4.9096e-06 | 
| "Pasteurellaceae"     | 3.9372   | 0.81793 | 1.4827e-06 | 1.4827e-05 | 
| "Prevotellaceae"      | 3.9769   | 0.84757 | 2.7043e-06 | 2.1635e-05 | 
| "Bifidobacteriaceae"  | 2.87     | 0.63484 | 6.1595e-06 | 4.1063e-05 | 
| "Enterobacteriaceae"  | 3.0798   | 0.85801 | 0.00033131 | 0.0018932  | 
| "Fusobacteriaceae"    | 3.8349   | 1.0842  | 0.00040474 | 0.0020237  | 
| "Micrococcaceae"      | 2.2196   | 0.79    | 0.0049611  | 0.022049   | 

| "Genus"                  | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|--------------------------|----------|---------|------------|------------| 
| "Veillonella"            | 5.1543   | 0.72825 | 1.4657e-12 | 9.8199e-11 | 
| "Shigella"               | 4.3745   | 0.84122 | 1.9909e-07 | 6.5364e-06 | 
| "Streptococcus"          | 3.1351   | 0.61137 | 2.9268e-07 | 6.5364e-06 | 
| "Haemophilus"            | 3.9829   | 0.81921 | 1.1631e-06 | 1.9482e-05 | 
| "Megasphaera"            | 5.4152   | 1.1367  | 1.8977e-06 | 2.5429e-05 | 
| "Prevotella"             | 4.1747   | 0.88327 | 2.2858e-06 | 2.5525e-05 | 
| "Fusobacterium"          | 4.1856   | 1.0564  | 7.4267e-05 | 0.00071084 | 
| "Fusicatenibacter"       | -2.3343  | 0.62475 | 0.00018676 | 0.0015641  | 
| "Rothia"                 | 2.3558   | 0.74219 | 0.001503   | 0.011189   | 
| "Candidatus_Soleaferrea" | -2.4742  | 0.8231  | 0.0026477  | 0.01774    | 

| "Species"                                                | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|----------------------------------------------------------|----------|---------|------------|------------| 
| "OTU1028-NN=Ruminococcus_gnavus_JN713312-D=95_6"         | 6.5111   | 0.93552 | 3.4073e-12 | 3.6976e-10 | 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"           | 5.1259   | 0.73858 | 3.9128e-12 | 3.6976e-10 | 
| "Streptococcus_infantis"                                 | 3.9275   | 0.62695 | 3.7421e-10 | 2.3575e-08 | 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89"    | -5.9611  | 0.99089 | 1.7888e-09 | 8.4519e-08 | 
| "Dialister_pneumosintes"                                 | 7.2606   | 1.2199  | 2.6535e-09 | 1.003e-07  | 
| "Streptococcus_anginosus"                                | 5.0796   | 0.86981 | 5.225e-09  | 1.6459e-07 | 
| "OTU790-NN=Barnesiella_intestinihominis_AB370251-D=84_6" | 6.7489   | 1.2191  | 3.093e-08  | 8.3511e-07 | 
| "OTU506-NN=Christensenella_minuta_AB490809-D=83_8"       | -5.7248  | 1.0396  | 3.6604e-08 | 8.6476e-07 | 
| "Veillonella_rogosae"                                    | 5.1187   | 0.95835 | 9.2351e-08 | 1.9394e-06 | 
| "Shigella_sonnei"                                        | 4.5891   | 0.87351 | 1.4918e-07 | 2.8195e-06 | 
| "OTU175-NN=Roseburia_hominis_AB661434-D=94_3"            | -4.8133  | 0.94715 | 3.7378e-07 | 6.4222e-06 | 
| "Parabacteroides_goldsteinii"                            | -4.698   | 0.97152 | 1.3269e-06 | 2.0898e-05 | 
| "Haemophilus_parainfluenzae"                             | 3.8547   | 0.85371 | 6.3254e-06 | 9.1961e-05 | 
| "OTU1339-NN=Clostridium_clariflavum_NR_102987_1-D=76_8"  | -4.704   | 1.1527  | 4.4876e-05 | 0.00060583 | 
| "OTU869-NN=Papillibacter_cinnamivorans_AF167711-D=90_2"  | -4.2336  | 1.1193  | 0.00015527 | 0.0019564  | 
| "OTU84-NN=Robinsoniella_peoriensis_AF445258-D=93_1"      | -2.952   | 0.7875  | 0.0001779  | 0.0021014  | 
| "OTU1000-NN=Bacteroides_vulgatus_BNRRR16SA-D=96_9"       | 2.3596   | 0.63638 | 0.00020908 | 0.0023244  | 
| "OTU1074-NN=Papillibacter_cinnamivorans_AF167711-D=85_3" | -4.5432  | 1.238   | 0.00024261 | 0.0025474  | 
| "OTU1678-NN=Veillonella_dispar_GQ422726-D=94_9"          | 5.4823   | 1.5065  | 0.00027357 | 0.0026489  | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6"             | 3.5841   | 0.98659 | 0.00028031 | 0.0026489  | 
| "Phascolarctobacterium_succinatutens"                    | 4.3877   | 1.2269  | 0.00034866 | 0.0031379  | 
| "OTU1033-NN=Ruminococcus_bromii_DQ882649-D=88"           | -2.8326  | 0.82496 | 0.00059558 | 0.0051166  | 
| "OTU1683-NN=Ruminococcus_bromii_DQ882649-D=90"           | -3.5819  | 1.0514  | 0.00065738 | 0.005402   | 
| "Fusicatenibacter_saccharivorans"                        | -2.2513  | 0.67178 | 0.00080451 | 0.0063356  | 
| "OTU1240-NN=OscillospiraClostridium_viride_X81125-D=86"  | -3.9442  | 1.1912  | 0.00092952 | 0.0070272  | 
| "OTU555-NN=Christensenella_minuta_AB490809-D=88_7"       | -3.7926  | 1.1618  | 0.0010972  | 0.0079761  | 
| "OTU867-NN=Christensenella_minuta_AB490809-D=91_9"       | -2.6129  | 0.83408 | 0.0017323  | 0.012126   | 
| "OTU1524-NN=Anaerofilum_pentosovorans_X97852-D=90_6"     | -2.8035  | 0.89844 | 0.0018058  | 0.012189   | 
| "Ruminococcus_callidus"                                  | -3.2155  | 1.0632  | 0.0024917  | 0.016239   | 
| "OTU855-NN=Roseburia_faecis_AY804149-D=95_5"             | -1.8723  | 0.63147 | 0.0030269  | 0.01907    | 
| "OTU252-NN=Rothia_aeria_AB753461-D=96_6"                 | 2.4169   | 0.81861 | 0.0031532  | 0.019224   | 
| "Granulicatella_adiacens"                                | 1.8893   | 0.64715 | 0.0035068  | 0.020712   | 
| "OTU342-NN=Clostridium_clariflavum_NR_102987_1-D=77_3"   | -2.4839  | 0.88497 | 0.0050033  | 0.028655   | 
| "Not_Assigned"                                           | 0.98141  | 0.35954 | 0.0063406  | 0.035246   | 
| "OTU962-NN=Oscillibacter_valericigenes_AB238598-D=89_7"  | -3.4972  | 1.2979  | 0.0070516  | 0.038079   | 
| "OTU379-NN=Dialister_succinatiphilus_AB370249-D=96_3"    | -3.4533  | 1.2983  | 0.0078173  | 0.041041   | 
| "Eggerthella_lenta"                                      | 2.1472   | 0.8134  | 0.0082967  | 0.042381   | 
| "OTU661-NN=Blautia_wexlerae_EF036467-D=94_9"             | -1.5468  | 0.58839 | 0.0085675  | 0.042612   | 
| "OTU1689-NN=Soleaferrea_massiliensis_JX101688-D=88_9"    | -2.3464  | 0.89842 | 0.0090102  | 0.043665   | 
| "OTU197-NN=Ruminococcus_albus_AY445596-D=93_7"           | -2.5669  | 0.99553 | 0.0099252  | 0.046896   | 

| "OTU"                                                    | "log2FC" | "lfcSE" | "Pvalues"  | "FDR"      | 
|----------------------------------------------------------|----------|---------|------------|------------| 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"           | 6.2078   | 0.77191 | 8.8328e-16 | 2.641e-13  | 
| "Streptococcus_1"                                        | 5.5123   | 0.69674 | 2.5426e-15 | 3.8012e-13 | 
| "OTU1028-NN=Ruminococcus_gnavus_JN713312-D=95_6"         | 6.4167   | 0.92343 | 3.6845e-12 | 3.6722e-10 | 
| "Bacteroides_3"                                          | -6.7879  | 1.0228  | 3.2156e-11 | 2.4036e-09 | 
| "Streptococcus_infantis"                                 | 3.8191   | 0.60326 | 2.4411e-10 | 1.4598e-08 | 
| "Dialister_pneumosintes"                                 | 7.2558   | 1.2099  | 2.0108e-09 | 1.0021e-07 | 
| "OTU506-NN=Christensenella_minuta_AB490809-D=83_8"       | -6.0393  | 1.0339  | 5.1736e-09 | 2.2099e-07 | 
| "Streptococcus_anginosus"                                | 4.9493   | 0.86599 | 1.0958e-08 | 4.0957e-07 | 
| "OTU790-NN=Barnesiella_intestinihominis_AB370251-D=84_6" | 6.7883   | 1.2109  | 2.0705e-08 | 6.2652e-07 | 
| "Streptococcus"                                          | 3.5874   | 0.64015 | 2.0954e-08 | 6.2652e-07 | 
| "Bacteroidales_5"                                        | -5.8559  | 1.0889  | 7.5456e-08 | 1.9139e-06 | 
| "Shigella_sonnei"                                        | 4.6813   | 0.87102 | 7.6814e-08 | 1.9139e-06 | 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89_2"  | -5.6423  | 1.0634  | 1.1205e-07 | 2.5772e-06 | 
| "OTU175-NN=Roseburia_hominis_AB661434-D=94_3"            | -4.8783  | 0.93763 | 1.9637e-07 | 4.1938e-06 | 
| "Paraprevotella_clara"                                   | 5.563    | 1.0745  | 2.2498e-07 | 4.4847e-06 | 
| "Parabacteroides_goldsteinii"                            | -4.7466  | 0.96598 | 8.9367e-07 | 1.67e-05   | 
| "Megasphaera"                                            | 5.4131   | 1.1279  | 1.593e-06  | 2.7819e-05 | 
| "Veillonella_rogosae"                                    | 4.3818   | 0.91495 | 1.6747e-06 | 2.7819e-05 | 
| "Clostridiales_23"                                       | -5.5476  | 1.1958  | 3.4947e-06 | 5.4995e-05 | 
| "Haemophilus_parainfluenzae"                             | 3.8578   | 0.8441  | 4.8698e-06 | 7.2804e-05 | 
| "Lachnospiraceae_7"                                      | -4.4495  | 0.9944  | 7.6574e-06 | 0.00010903 | 
| "Bacteroides_1"                                          | 5.2738   | 1.209   | 1.288e-05  | 0.00017506 | 
| "OTU1000-NN=Bacteroides_vulgatus_BNRRR16SA-D=96_9"       | 2.6943   | 0.62703 | 1.7315e-05 | 0.00022509 | 
| "OTU1339-NN=Clostridium_clariflavum_NR_102987_1-D=76_8"  | -4.8767  | 1.142   | 1.9523e-05 | 0.00024322 | 
| "Lachnospiraceae"                                        | -2.2711  | 0.54315 | 2.898e-05  | 0.0003466  | 
| "OTU84-NN=Robinsoniella_peoriensis_AF445258-D=93_1"      | -3.0892  | 0.76335 | 5.1902e-05 | 0.00059688 | 
| "Pasteurellaceae"                                        | 3.6134   | 0.91724 | 8.1674e-05 | 0.00090447 | 
| "unclassified_39"                                        | -3.8967  | 1.0098  | 0.0001139  | 0.0012162  | 
| "Fusobacterium"                                          | 4.1157   | 1.0745  | 0.00012799 | 0.0013196  | 
| "Oscillospira_4"                                         | -4.1742  | 1.1007  | 0.0001493  | 0.0014729  | 
| "Prevotella_1"                                           | 4.0681   | 1.0743  | 0.00015271 | 0.0014729  | 
| "Bacteria_5"                                             | -3.7444  | 1.0125  | 0.00021709 | 0.0020039  | 
| "OTU1074-NN=Papillibacter_cinnamivorans_AF167711-D=85_3" | -4.5908  | 1.2429  | 0.00022117 | 0.0020039  | 
| "Clostridiales_54"                                       | -3.6997  | 1.0256  | 0.00030932 | 0.0026638  | 
| "Fusicatenibacter_saccharivorans"                        | -2.3442  | 0.65022 | 0.00031182 | 0.0026638  | 
| "OTU869-NN=Papillibacter_cinnamivorans_AF167711-D=90_2"  | -3.9957  | 1.1141  | 0.00033526 | 0.0027836  | 
| "unclassified_5"                                         | -4.0589  | 1.134   | 0.00034446 | 0.0027836  | 
| "OTU1678-NN=Veillonella_dispar_GQ422726-D=94_9"          | 5.2826   | 1.4926  | 0.00040142 | 0.0031586  | 
| "unclassified_4"                                         | -3.2215  | 0.92199 | 0.00047579 | 0.0036477  | 
| "Phascolarctobacterium_succinatutens"                    | 4.1989   | 1.2227  | 0.0005943  | 0.0044424  | 
| "OTU516-NN=Prevotella_copri_AB064923-D=96_6"             | 3.4175   | 0.99732 | 0.00061111 | 0.0044567  | 
| "OTU1240-NN=OscillospiraClostridium_viride_X81125-D=86"  | -4.0036  | 1.1974  | 0.00082714 | 0.0058885  | 
| "Clostridiales_20"                                       | -3.051   | 0.91538 | 0.00085911 | 0.0059738  | 
| "OTU197-NN=Ruminococcus_albus_AY445596-D=93_7"           | -3.3576  | 1.0119  | 0.00090675 | 0.0061618  | 
| "OTU1524-NN=Anaerofilum_pentosovorans_X97852-D=90_6"     | -2.9303  | 0.89419 | 0.0010488  | 0.0068981  | 
| "OTU1683-NN=Ruminococcus_bromii_DQ882649-D=90"           | -3.4039  | 1.0397  | 0.0010613  | 0.0068981  | 
| "OTU555-NN=Christensenella_minuta_AB490809-D=88_7_1"     | -3.7906  | 1.1655  | 0.0011441  | 0.0072787  | 
| "Barnesiella_10"                                         | -3.5754  | 1.1089  | 0.0012624  | 0.0078635  | 
| "Ruminococcus_callidus"                                  | -3.3784  | 1.0654  | 0.0015186  | 0.0090836  | 
| "unclassified_38"                                        | -3.0702  | 0.96819 | 0.001519   | 0.0090836  | 
| "OTU867-NN=Christensenella_minuta_AB490809-D=91_9"       | -2.595   | 0.82302 | 0.0016162  | 0.0093557  | 
| "Paraprevotella_clara_1"                                 | 4.1391   | 1.3136  | 0.0016271  | 0.0093557  | 
| "OTU252-NN=Rothia_aeria_AB753461-D=96_6"                 | 2.4406   | 0.7932  | 0.0020918  | 0.011801   | 
| "Bacteria_3"                                             | 1.0919   | 0.35854 | 0.0023233  | 0.012729   | 
| "Bacteria"                                               | 1.0506   | 0.34523 | 0.0023415  | 0.012729   | 
| "Bacteria_8"                                             | 1.0377   | 0.34403 | 0.0025595  | 0.013666   | 
| "OTU855-NN=Roseburia_faecis_AY804149-D=95_5"             | -1.8659  | 0.62198 | 0.0027009  | 0.014168   | 
| "OTU342-NN=Clostridium_clariflavum_NR_102987_1-D=77_3"   | -2.6048  | 0.8784  | 0.003023   | 0.015584   | 
| "Bacteria_6"                                             | 1.0087   | 0.34363 | 0.0033307  | 0.016879   | 
| "Barnesiella_intestinihominis_1"                         | -2.6278  | 0.90084 | 0.0035341  | 0.017559   | 
| "Lachnoclostridium_1"                                    | -2.04    | 0.70035 | 0.0035823  | 0.017559   | 
| "Granulicatella_adiacens"                                | 1.7855   | 0.63206 | 0.0047308  | 0.022815   | 
| "OTU962-NN=Oscillibacter_valericigenes_AB238598-D=89_7"  | -3.6398  | 1.2994  | 0.0050939  | 0.024176   | 
| "Eubacterium_2"                                          | 1.8744   | 0.68114 | 0.0059252  | 0.02734    | 
| "unclassified_3"                                         | -2.2028  | 0.80077 | 0.0059434  | 0.02734    | 
| "Bacteroides"                                            | 1.8081   | 0.67118 | 0.0070634  | 0.031656   | 
| "Atopobium"                                              | 2.0837   | 0.77393 | 0.0070936  | 0.031656   | 
| "OTU379-NN=Dialister_succinatiphilus_AB370249-D=96_3"    | -3.4238  | 1.3097  | 0.0089411  | 0.039315   | 
| "Erysipelatoclostridium"                                 | 2.1258   | 0.83996 | 0.011379   | 0.049308   | 
| "OTU661-NN=Blautia_wexlerae_EF036467-D=94_9"             | -1.4804  | 0.58667 | 0.011624   | 0.049649   | 

### Examples:

 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Bifidobacteriaceae_1.png)

 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Desulfovibrionaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Enterobacteriaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Fusobacteriaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Micrococcaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Pasteurellaceae_1.png)

 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Prevotellaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Streptococcaceae_1.png)
 
 ![alt text](/MicrobiomAnalysis/Liver-all/DESeq2/Veillonellaceae_1.png)

## LEfSe

The LEfSe algorithm was developed for biomarker discovery and interpretation for metagenomics data. It employs Kruskal- Wallis rank sum test to detect features with significant differential abundance with regard to class labels, followed by Linear Discriminant Analysis to evaluate the relevance or effect size of differential abundant features.

See results in folder:

 ![LEfSe-Results](/MicrobiomAnalysis/Liver-all/LEfSe/)

| "Class"               | "Pvalues" | "FDR"    | "healthy"       | "normal"        | "fibrose"       | "zirrhose"      | "LDAscore" | 
|-----------------------|-----------|----------|------------------|------------------|------------------|------------------|------------| 
| "Deltaproteobacteria" | 0.0013122 | 0.018371 | 14155.9267153139 | 5156.86054950193 | 947.58814229856  | 799.295779140957 | 3.82       | 
| "Verrucomicrobiae"    | 0.0067136 | 0.046995 | 181624.740121624 | 87312.3305394486 | 1501472.27110865 | 87681.8900794324 | 5.85       | 

| "Order"              | "Pvalues" | "FDR"    | "healthy"       | "normal"        | "fibrose"      | "zirrhose"      | "LDAscore" | 
|----------------------|-----------|----------|------------------|------------------|-----------------|------------------|------------| 
| "Desulfovibrionales" | 0.0013122 | 0.030181 | 14155.9267153139 | 5156.86054950193 | 947.58814229856 | 799.295779140957 | 3.82       | 

| "Species"                                                    | "Pvalues"  | "FDR"     | "healthy"       | "normal"        | "fibrose"       | "zirrhose"      | "LDAscore" | 
|--------------------------------------------------------------|------------|-----------|------------------|------------------|------------------|------------------|------------| 
| "OTU1337-NN=Desulfovibrio_desulfuricans_DVURRDA-D=89"        | 2.5068e-05 | 0.0047628 | 8462.44760624677 | 2510.48286950075 | 488.494334667348 | 43.9942737326595 | 3.62       | 
| "OTU506-NN=Christensenella_minuta_AB490809-D=83_8"           | 6.2399e-05 | 0.0048236 | 97187.570906333  | 52045.1492756922 | 7286.47990258686 | 475.00243296263  | 4.68       | 
| "OTU916-NN=Eubacterium_siraeum_EUBRRDV-D=93_2"               | 7.6162e-05 | 0.0048236 | 17393.7003925028 | 2761.75765933834 | 21.4443636706746 | 3868.71050294625 | 3.94       | 
| "OTU383-NN=Clostridium_populeti_X71853-D=91_5"               | 0.00023753 | 0.011283  | 14734.2585420822 | 35877.5747128127 | 3041.57185887673 | 9023.65047746983 | 4.22       | 
| "OTU1033-NN=Ruminococcus_bromii_DQ882649-D=88"               | 0.0003066  | 0.011651  | 2375.2152708857  | 260.373731840983 | 1194.41725549619 | 711.329499727706 | 3.02       | 
| "OTU1335-NN=Oscillibacter_valericigenes_AB238598-D=95_6"     | 0.00043939 | 0.013914  | 3898.11656376373 | 3840.73352562476 | 1595.24954199799 | 1483.65348546604 | 3.08       | 
| "Dialister_pneumosintes"                                     | 0.00060686 | 0.015669  | 0                | 2195.80314308171 | 285.697562754108 | 43077.0779944808 | 4.33       | 
| "OTU661-NN=Blautia_wexlerae_EF036467-D=94_9"                 | 0.00065975 | 0.015669  | 28024.8817821196 | 15003.9107284982 | 3482.76886188742 | 4462.72756529591 | 4.09       | 
| "OTU694-NN=Clostridium_clariflavum_NR_102987_1-D=85_2"       | 0.0010549  | 0.02227   | 9440.38516187504 | 6069.24180458566 | 3966.53325813013 | 2198.76804721345 | 3.56       | 
| "OTU293-NN=Clostridium_glycyrrhizinilyticum_AB233029-D=91_2" | 0.0011866  | 0.022546  | 9076.64778195034 | 3506.21944381522 | 30.7795070376657 | 2307.27887656725 | 3.66       | 
| "OTU1689-NN=Soleaferrea_massiliensis_JX101688-D=88_9"        | 0.0013189  | 0.022782  | 5231.73992941993 | 3051.8222455415  | 3475.82665060504 | 560.277041563892 | 3.37       | 
| "OTU342-NN=Clostridium_clariflavum_NR_102987_1-D=77_3"       | 0.0015065  | 0.023852  | 17644.4159246064 | 8520.08840776765 | 3678.86865878213 | 1450.10859534348 | 3.91       | 
| "OTU1482-NN=Veillonella_atypica_X84007-D=96_9"               | 0.0023425  | 0.033569  | 114489.215468368 | 244072.179110858 | 237031.899462028 | 1410885.47779397 | 5.81       | 
| "OTU1339-NN=Clostridium_clariflavum_NR_102987_1-D=76_8"      | 0.0024735  | 0.033569  | 13615.7631275547 | 4655.63942077853 | 30733.8974395363 | 73.8525195105683 | 4.19       | 
| "OTU1240-NN=OscillospiraClostridium_viride_X81125-D=86"      | 0.0029348  | 0.036122  | 1423.02412528051 | 112.232582722412 | 747.994393800127 | 0                | 2.85       | 
| "Parabacteroides_goldsteinii"                                | 0.0030418  | 0.036122  | 11720.6757435959 | 2750.28181220383 | 3461.66709320441 | 209.726669436687 | 3.76       | 
| "OTU902-NN=Prevotella_loescheii_AB547688-D=88_8"             | 0.0039699  | 0.04437   | 8848.44554766887 | 1727.23649559518 | 4682.19080552405 | 1333.62736346102 | 3.58       | 

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-all/LEfSe/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-all/LEfSe/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-all/LEfSe/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-all/LEfSe/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-all/LEfSe/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-all/LEfSe/7_1.png)


## Random Forests

Random Forests is a powerful machine learning algorithms for classification as well as for identification of predictive features (biomakers). It operates by constructing a multitude of decision trees ("forests") at training time and predicting the class as the majority vote of the individual trees.

### Phylum

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c2_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c3_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c4_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c5_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c6_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c7_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f7_1.png)

### OTU

 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/c8_1.png)
 ![alt text](/MicrobiomAnalysis/Liver-all/RandomForests/f8_1.png)

## Core Microbiome

Core microbiome refers to the set of taxa that are detected in a high fraction of the population above a given abundance threshold. The count data is transformed to compositional (relative) abundances in order to perform such analysis.

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-all/CoreMicrobiome/genus_1.png)

## Hierarchical Clustering & Heatmap Visualization

### Genus

 ![alt text](/MicrobiomAnalysis/Liver-all/Heatmap/genus_2.png)
 
## Dendrogram (OTUs)

 ![alt text](/MicrobiomAnalysis/Liver-all/Dendrogram/OTU_1.png)

## Correlation Analysis (Genus)

 ![alt text](/MicrobiomAnalysis/Liver-all/CorrelationAnalysis/genus_2.png)

