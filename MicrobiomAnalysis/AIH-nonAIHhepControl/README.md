# AIH versus nonAIHhepControl

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/norm_libsizes_0_1.png)

## Taxonomic Bar Plots

by Sample:

### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/7_1.png)

by Group

### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/g2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/g3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/TaxBars/g4_1.png)


## Alpha Diversity

Observed	p-value: 0.0078032; [T-test] statistic: 2.8509

Observed	p-value: 0.0097489; [Mann-Whitney] statistic: 275

Chao1	p-value: 0.0075928; [T-test] statistic: 2.8695

Chao1	p-value: 0.0053966; [Mann-Whitney] statistic: 280

ACE	p-value: 0.0045493; [T-test] statistic: 3.0659

ACE	p-value: 0.0036375; [Mann-Whitney] statistic: 284

Shannon	p-value: 0.04944; [T-test] statistic: 2.037

Shannon	p-value: 0.056954; [Mann-Whitney] statistic: 251

Simpson	p-value: 0.091489; [T-test] statistic: 1.7327

Simpson	p-value: 0.084115; [Mann-Whitney] statistic: 245

Fisher	p-value: 0.006944; [T-test] statistic: 2.9136

Fisher	p-value: 0.0036375; [Mann-Whitney] statistic: 284

### Observed

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/Observed_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/gObserved_1.png)

### Chao1

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/Chao1_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/gChao1_1.png)

### ACE

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/ACE_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/gACE_1.png)

### Shannon

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/Shannon_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/gShannon_1.png)

### Simpson

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/Simpson_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/gSimpson_1.png)

### Fisher

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/Fisher_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Alphadiversity/gFisher_1.png)


## Beta Diversity

The statistical methods test the strength and statistical significance of sample groupings based on distance matrix. Users can choose from different statistical methods such as:

Permutational multivariate analysis of variance (PERMANOVA), is a non-parametric multivariate statistical test. PERMANOVA is used to compare groups of objects and test the null hypothesis that the centroids and dispersion of the groups as defined by measure space are equivalent for all groups.

adonis/ANOSIM: Tests whether two or more groups of samples are significantly different based on a categorical variable found in the sample mapping file. An R value near 1 means that there is dissimilarity between the groups, while an R value near 0 indicates no significant dissimilarity between the groups.

PERMDISP: This method analyzes the multivariate homogeneity of group dispersions (variances).This test is different from two others in that it specifically tests for differences in the spread (dispersion, variability) among groups. In essence, it determines whether the variances of groups of samples are significantly different.

[PERMANOVA] R-squared: 0.038991; p-value < 0.063

[ANOSIM] R: -0.025411; p-value < 0.685

[PERMDISP] F-value: 7.4332; p-value: 0.0097279

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Betadiversity/PCoA_1.png)
 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/Betadiversity/MNDS_1.png)

## Univariate Analysis
 
| "Phylum"          | "Pvalues" | "FDR"   | "Statistics" | 
|-------------------|-----------|---------|--------------| 
| "Bacteroidetes"   | 0.034735  | 0.24315 | 258          | 
| "Firmicutes"      | 0.16805   | 0.58816 | 135          | 
| "Proteobacteria"  | 0.26283   | 0.61326 | 144          | 
| "Actinobacteria"  | 0.40408   | 0.68431 | 214          | 
| "Unclassified"    | 0.57399   | 0.68431 | 203.5        | 
| "Fusobacteria"    | 0.58655   | 0.68431 | 166          | 
| "Verrucomicrobia" | 0.8304    | 0.8304  | 192          | 

| "Class"                 | "Pvalues" | "FDR"   | "Statistics" | 
|-------------------------|-----------|---------|--------------| 
| "Negativicutes"         | 0.020251  | 0.24835 | 103          | 
| "Bacteroidia"           | 0.034735  | 0.24835 | 258          | 
| "Erysipelotrichia"      | 0.04967   | 0.24835 | 253          | 
| "Epsilonproteobacteria" | 0.095417  | 0.33778 | 134.5        | 
| "Gammaproteobacteria"   | 0.12059   | 0.33778 | 129          | 
| "Clostridia"            | 0.13511   | 0.33778 | 237          | 
| "Actinobacteria"        | 0.40408   | 0.67679 | 214          | 
| "Deltaproteobacteria"   | 0.52158   | 0.67679 | 206          | 
| "Betaproteobacteria"    | 0.52502   | 0.67679 | 207          | 
| "Coriobacteriia"        | 0.52907   | 0.67679 | 202          | 
| "Alphaproteobacteria"   | 0.54366   | 0.67679 | 206          | 
| "unclassified"          | 0.57399   | 0.67679 | 203.5        | 
| "Fusobacteriia"         | 0.58655   | 0.67679 | 166          | 
| "Verrucomicrobiae"      | 0.8304    | 0.85452 | 192          | 
| "Bacilli"               | 0.85452   | 0.85452 | 191          | 

| "Order"                   | "Pvalues" | "FDR"   | "Statistics" | 
|----------------------|-----------|---------|--------------| 
| "Selenomonadales"    | 0.020251  | 0.36425 | 103          | 
| "Bacteroidales"      | 0.034735  | 0.36425 | 258          | 
| "Erysipelotrichales" | 0.04967   | 0.36425 | 253          | 
| "Campylobacterales"  | 0.095417  | 0.46909 | 134.5        | 
| "Enterobacteriales"  | 0.11961   | 0.46909 | 129          | 
| "Clostridiales"      | 0.13511   | 0.46909 | 237          | 
| "Actinomycetales"    | 0.14925   | 0.46909 | 235          | 
| "Caulobacterales"    | 0.29732   | 0.81763 | 221          | 
| "Pasteurellales"     | 0.46612   | 0.85452 | 210          | 
| "Desulfovibrionales" | 0.52158   | 0.85452 | 206          | 
| "Burkholderiales"    | 0.52502   | 0.85452 | 207          | 
| "Pseudomonadales"    | 0.53815   | 0.85452 | 206          | 
| "unclassified"       | 0.57399   | 0.85452 | 203.5        | 
| "Sphingomonadales"   | 0.58187   | 0.85452 | 204          | 
| "Fusobacteriales"    | 0.58655   | 0.85452 | 166          | 
| "Bifidobacteriales"  | 0.68237   | 0.85452 | 199          | 
| "Bacillales"         | 0.70123   | 0.85452 | 197          | 
| "Gemellales"         | 0.73174   | 0.85452 | 196          | 
| "Coriobacteriales"   | 0.81049   | 0.85452 | 193          | 
| "Rhizobiales"        | 0.8199    | 0.85452 | 192          | 
| "Verrucomicrobiales" | 0.8304    | 0.85452 | 192          | 
| "Lactobacillales"    | 0.85452   | 0.85452 | 191          | 

| "Family"                | "Pvalues" | "FDR"   | "Statistics" | 
|-------------------------|-----------|---------|--------------| 
| "Oscillospiraceae"      | 0.011252  | 0.43885 | 272          | 
| "Peptostreptococcaceae" | 0.025665  | 0.44424 | 262          | 
| "Veillonellaceae"       | 0.043166  | 0.44424 | 113          | 
| "Erysipelotrichaceae"   | 0.04967   | 0.44424 | 253          | 
| "Rikenellaceae"         | 0.056954  | 0.44424 | 251          | 
| "Prevotellaceae"        | 0.074112  | 0.44795 | 247          | 
| "Campylobacteraceae"    | 0.095417  | 0.44795 | 134.5        | 
| "Actinomycetaceae"      | 0.11268   | 0.44795 | 240          | 
| "Enterobacteriaceae"    | 0.11961   | 0.44795 | 129          | 
| "Porphyromonadaceae"    | 0.1277    | 0.44795 | 238          | 
| "unclassified"          | 0.1277    | 0.44795 | 238          | 
| "Mogibacteriaceae"      | 0.14178   | 0.44795 | 235          | 
| "Lactobacillaceae"      | 0.14932   | 0.44795 | 235          | 
| "Ruminococcaceae"       | 0.26283   | 0.68808 | 224          | 
| "Bacteroidaceae"        | 0.28801   | 0.68808 | 222          | 
| "Caulobacteraceae"      | 0.29732   | 0.68808 | 221          | 
| "Christensenellaceae"   | 0.32053   | 0.68808 | 219          | 
| "Acidaminococcaceae"    | 0.33081   | 0.68808 | 218.5        | 
| "Clostridiaceae"        | 0.34299   | 0.68808 | 218          | 
| "Lachnospiraceae"       | 0.35769   | 0.68808 | 217          | 
| "Moraxellaceae"         | 0.3705    | 0.68808 | 214          | 
| "Sutterellaceae"        | 0.43104   | 0.76411 | 212          | 
| "Pasteurellaceae"       | 0.46612   | 0.79037 | 210          | 
| "Desulfovibrionaceae"   | 0.52158   | 0.84757 | 206          | 
| "Sphingomonadaceae"     | 0.58187   | 0.87529 | 204          | 
| "Fusobacteriaceae"      | 0.58655   | 0.87529 | 166          | 
| "Oxalobacteraceae"      | 0.63722   | 0.87529 | 201          | 
| "Micrococcaceae"        | 0.66992   | 0.87529 | 169          | 
| "Pseudomonadaceae"      | 0.67701   | 0.87529 | 199          | 
| "Bifidobacteriaceae"    | 0.68237   | 0.87529 | 199          | 
| "Staphylococcaceae"     | 0.70123   | 0.87529 | 197          | 
| "Gemellaceae"           | 0.73174   | 0.87529 | 196          | 
| "Eubacteriaceae"        | 0.74557   | 0.87529 | 196          | 
| "Comamonadaceae"        | 0.80821   | 0.87529 | 193          | 
| "Coriobacteriaceae"     | 0.81049   | 0.87529 | 193          | 
| "Bradyrhizobiaceae"     | 0.8199    | 0.87529 | 192          | 
| "Verrucomicrobiaceae"   | 0.8304    | 0.87529 | 192          | 
| "Carnobacteriaceae"     | 0.92038   | 0.9446  | 188          | 
| "Streptococcaceae"      | 0.96625   | 0.96625 | 186          | 

| "Geuns"                  | "Pvalues" | "FDR"   | "Statistics" | 
|--------------------------|-----------|---------|--------------| 
| "Aggregatibacter"        | 0.010702  | 0.17487 | 259.5        | 
| "Turicibacter"           | 0.01116   | 0.17487 | 273          | 
| "Blautia"                | 0.011252  | 0.17487 | 272          | 
| "Oscillibacter"          | 0.011252  | 0.17487 | 272          | 
| "Dorea"                  | 0.013426  | 0.17487 | 271          | 
| "Terrisporobacter"       | 0.015206  | 0.17487 | 269          | 
| "Barnesiella"            | 0.026882  | 0.26498 | 262          | 
| "Intestinibacter"        | 0.041163  | 0.3071  | 256          | 
| "Odoribacter"            | 0.043749  | 0.3071  | 255          | 
| "Alistipes"              | 0.056954  | 0.3071  | 251          | 
| "Fusicatenibacter"       | 0.057476  | 0.3071  | 251          | 
| "Butyricimonas"          | 0.057719  | 0.3071  | 250          | 
| "Faecalibacterium"       | 0.06508   | 0.3071  | 249          | 
| "Eggerthella"            | 0.065541  | 0.3071  | 123.5        | 
| "Lactococcus"            | 0.06676   | 0.3071  | 246          | 
| "Prevotella"             | 0.095154  | 0.36855 | 243          | 
| "Campylobacter"          | 0.095417  | 0.36855 | 134.5        | 
| "Parabacteroides"        | 0.10608   | 0.36855 | 241          | 
| "Shigella"               | 0.10653   | 0.36855 | 127          | 
| "Paraprevotella"         | 0.11075   | 0.36855 | 238          | 
| "Actinomyces"            | 0.11268   | 0.36855 | 240          | 
| "Brevundimonas"          | 0.11751   | 0.36855 | 235.5        | 
| "Anaerofilum"            | 0.13879   | 0.41119 | 228          | 
| "Parasutterella"         | 0.14302   | 0.41119 | 234          | 
| "Lactobacillus"          | 0.14932   | 0.41211 | 235          | 
| "Gemmiger"               | 0.1724    | 0.45752 | 232          | 
| "Intestinimonas"         | 0.18666   | 0.47701 | 228          | 
| "Veillonella"            | 0.27523   | 0.67612 | 145          | 
| "Bacteroides"            | 0.28801   | 0.67612 | 222          | 
| "Christensenella"        | 0.29691   | 0.67612 | 220.5        | 
| "Coprococcus"            | 0.30376   | 0.67612 | 220.5        | 
| "unclassified"           | 0.31473   | 0.67617 | 220          | 
| "Phascolarctobacterium"  | 0.33081   | 0.67617 | 218.5        | 
| "Kluyvera"               | 0.33364   | 0.67617 | 152          | 
| "Clostridium"            | 0.34299   | 0.67617 | 218          | 
| "Acinetobacter"          | 0.3705    | 0.70582 | 214          | 
| "Candidatus_Soleaferrea" | 0.37848   | 0.70582 | 215          | 
| "Phenylobacterium"       | 0.3996    | 0.71039 | 214          | 
| "Roseburia"              | 0.40408   | 0.71039 | 214          | 
| "Slackia"                | 0.41182   | 0.71039 | 208          | 
| "Sutterella"             | 0.42715   | 0.71887 | 212          | 
| "Collinsella"            | 0.48298   | 0.79347 | 209          | 
| "Haemophilus"            | 0.5018    | 0.80522 | 208          | 
| "Bilophila"              | 0.52158   | 0.81124 | 206          | 
| "Senegalemassilia"       | 0.52907   | 0.81124 | 202          | 
| "Sphingomonas"           | 0.58187   | 0.84077 | 204          | 
| "Fusobacterium"          | 0.58655   | 0.84077 | 166          | 
| "Allisonella"            | 0.62959   | 0.84077 | 199          | 
| "Erysipelatoclostridium" | 0.6312    | 0.84077 | 201          | 
| "Massilia"               | 0.63722   | 0.84077 | 201          | 
| "Ruminococcus"           | 0.64138   | 0.84077 | 201          | 
| "Anaerotruncus"          | 0.64813   | 0.84077 | 198.5        | 
| "Atopobium"              | 0.6679    | 0.84077 | 199          | 
| "Rothia"                 | 0.66992   | 0.84077 | 169          | 
| "Pseudomonas"            | 0.67701   | 0.84077 | 199          | 
| "Bifidobacterium"        | 0.68237   | 0.84077 | 199          | 
| "Staphylococcus"         | 0.70123   | 0.84886 | 197          | 
| "Gemella"                | 0.73174   | 0.87052 | 196          | 
| "Oscillospira"           | 0.76704   | 0.89528 | 195          | 
| "Lachnoclostridium"      | 0.78619   | 0.89528 | 194          | 
| "Pelomonas"              | 0.80821   | 0.89528 | 193          | 
| "Dialister"              | 0.81049   | 0.89528 | 175          | 
| "Bradyrhizobium"         | 0.8199    | 0.89528 | 192          | 
| "Akkermansia"            | 0.8304    | 0.89528 | 192          | 
| "Granulicatella"         | 0.92038   | 0.97702 | 188          | 
| "Megasphaera"            | 0.94225   | 0.98508 | 181          | 
| "Eubacterium"            | 0.96625   | 0.98875 | 186          | 
| "Tyzzerella"             | 0.97718   | 0.98875 | 182.5        | 
| "Streptococcus"          | 0.98875   | 0.98875 | 185          | 



## metagenomeSeq

none

## DDseq2

see folder

## LEfSe

### Phylum

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/LEfSe/2_1.png)

### Class

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/LEfSe/3_1.png)

### Order

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/LEfSe/4_1.png)

### Family

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/LEfSe/5_1.png)

### Genus

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/LEfSe/6_1.png)

### Species

 ![alt text](/MicrobiomAnalysis/AIH-nonAIHhepControl/LEfSe/7_1.png)

