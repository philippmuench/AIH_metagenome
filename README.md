# code and data to reproduce Münch _et al._ manuscript

Specific alterations of the intestinal microbiome in autoimmune hepatitis show partially liver disease-specific patterns being also related to the extent of liver parenchymatous tissue remodeling

# Table of Contents  
[Usage](#Usage)  
[Figure 1 (alpha diversity)](#figure-1)  
[Figure 2 (class level)](#figure-2)  
[Figure 3 (genus level)](#figure-3)  
[Figure 5 (alpha diveristy with hep)](#figure-5)  
[Figure 6 (marker correlation)](#figure-6)  
[Citation](#citation)

# Usage

clone this repo

```
git clone https://github.com/philippmuench/AIH_new.git
cd AIH_new
```
open R and install R packages please type these commands to your R console

```r
install.packages("reshape2")
install.packages("ggplot2")
install.packages("pander")
install.packages("ggtern")
install.packages("limma")
install.packages("biom")
install.packages("colorspace")
install.packages("grid")
install.packages("gridExtra")
install.packages("ggrepel")
source("https://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
biocLite("phyloseq")
```

# Figures
## Figure 1
**Description:**

> Figure 1: Mean and SD of alpha diversity of gut microbiota composition by disease status (AIH: black; non-AIH hepatitic diseases; grey; healthy: green)

**Notes:**

others group: NASH, alkohol, others, ALD, HCV

**File:**

- [figure 1 (observed OTUs index)](results/figure1/figure_1_index_observed.pdf)
- [figure 1 (chao1 index)](results/figure1/figure_1_index_chao1.pdf)
- [figure 1 (Shannon index)](results/figure1/figure_1_index_shannon.pdf)

**Statistics:**

- observed

```
    Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = value ~ type, data = df.adiv)

$type
                    diff        lwr        upr     p adj
control-AIH     143.8322   17.93869  269.72579 0.0215012
others-AIH     -131.6277 -252.41317  -10.84226 0.0296658
others-control -275.4600 -390.48441 -160.43550 0.0000011
```

- chao1
```
    Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = value ~ type, data = df.adiv)

$type
                    diff       lwr        upr     p adj
control-AIH     188.6858   14.3814  362.99019 0.0309829
others-AIH     -202.2373 -369.4693  -35.00522 0.0140803
others-control -390.9231 -550.1788 -231.66734 0.0000007
```

- shannon
```
  Using sample, type, type2 as id variables
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = value ~ type, data = df.adiv)

$type
                     diff       lwr        upr     p adj
control-AIH    -0.2994906 -1.252769 0.65378825 0.7308472
others-AIH     -0.8918631 -1.806463 0.02273673 0.0574292
others-control -0.5923725 -1.463349 0.27860441 0.2384276
```

**Reproduce:**

```r
source('figure1.R') # P values are written to stdout
```


## Figure 2

**Description:** 

> Figure2: a) Average relative abundance othe most relevant phyla (?) in health, AIH and control asrevealed by the 16S rRNA gene ribotyping. b) Ternary plot of all OTUs detected in the data set with RA > 1% in at least one sample. Each circle represent one class. The size of each circle represents its relative abundance (weighted average) of the OTUs associated with this class. The position of each circle is determined by the contribution of the indicated disease status to the total relative abundance. c) FDR adjusted P values of differently expressed bacterial classes between disese groups.

**Notes:**

Applied filtering: `removeRelAb()` with 0.01 (genera with rel. abundance < 1% are removed)

**File:**
- [figure 2 (merged)](results/figure2/figure_2.pdf)
- [figure 2a (bar plot)](results/figure2/figure_2_bar.pdf)
- [figure 2b up (triplot)](results/figure2/figure_2_tri.pdf)
- [figure 2b down (difference)](results/figure2/figure_2_diff.pdf)
- [data to reproduce figure 2c (ZIG analysis on class level, AIH vs. control)](results/figure2/figure_2_c_aih_vs_control.csv)
- [data to reproduce figure 2c (ZIG analysis on class level, AIH vs. healthy)](results/figure2/figure_2_c_aih_vs_healthy.csv)

**Reproduce:**

```r
source('figure2.R') # to reproduce figure 2 a, b
source('figure2_c.R') # to reproduce data shown in figure 2 c
```

## Figure 3

**Description:**

> Figure 3: b Ternary plot of all OTUs classified to class level. Each circle represents the mean abundance of OTUs associated to one class. The size of each circle represents its relative abundance (weighted average). The position of each circle is determined by the contribution of the indicated disease status to the total relative abundance. c) Heatmap of RA abundance of OTUs pooled to species level and single linkage hierarchical clustering.


**Files:**

- [data to reproduce figure 3a (ZIG analysis on genus level, AIH vs. control)](results/figure3/figure3_a_zig_aih_vs_control.tsv)
- [data to reproduce figure 3a (ZIG analysis on genus level, AIH vs. healthy)](results/figure3/figure3_a_zig_aih_vs_helathy.tsv)
- [data to reproduce figure 3a (ZIG analysis on genus level, AIH vs. healthy)](results/figure3/figure3_a_zig_healthy_vs_control.tsv)
- [figure 3 b](results/figure3/figure_3_b.pdf)
- [figure 3 c (heatmap)](results/figure3/figure_3_c_heatmap.pdf)
- [figure 3 c (dendogram row)](results/figure3/figure_3_c_dendogram_row.pdf)
- [figure 3 c (dendogram col)](results/figure3/figure_3_c_dendogram_col.pdf)

**Reproduce:**

```r
source('figure3_a.R') # generate tables with ZIG analysis results
source('figure3_b.R')
source('figure3_c.R')
```

## Figure 5

**Description:**

> Figure 5: Mean and SD of alpha diversity and hepatopathy status.

**File:**
- [figure 5](results/figure5/figure_5.pdf)


**Results:**
```
    Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = chao1 ~ marker, data = alpha)

$marker
                                                   diff        lwr        upr     p adj
hepatop_ohne_veraenderung-                    676.72081   20.24452 1333.19711 0.0402698
keine_hepatopathie-                           872.74112  249.17307 1496.30916 0.0020684
mit_fibrose-                                  608.94853  -25.85588 1243.75294 0.0659437
mit_zirrhose-                                 536.56865  -85.51294 1158.65024 0.1217344
keine_hepatopathie-hepatop_ohne_veraenderung  196.02030  -88.59827  480.63887 0.3070441
mit_fibrose-hepatop_ohne_veraenderung         -67.77229 -376.23179  240.68721 0.9711754
mit_zirrhose-hepatop_ohne_veraenderung       -140.15217 -421.49915  141.19482 0.6261394
mit_fibrose-keine_hepatopathie               -263.79259 -494.06035  -33.52483 0.0171538
mit_zirrhose-keine_hepatopathie              -336.17247 -528.60974 -143.73519 0.0000800
mit_zirrhose-mit_fibrose                      -72.37988 -298.59136  153.83161 0.8942774
```

**Reproduce:**

```r
source('figure5.R') # P values are written to stdout
```

## Figure 6

**Description:**

> Figure 6: Correlation between alpha diversity and abundance of indicated marker.

**File:**

- [figure 6](results/figure6/figure_6_all.pdf)
- [figure 6 (without all)](results/figure6/figure_6.pdf)

**Results:**

Correlation coefficient:

shannon

| Marker | group  | cor | P value  |
| --------- | ------------- |:-------------:| -----:|
| ifap | AIH | 0.2498182 | 0.3508 |
| ifap | healthy | -0.04517274 | 0.8891 |
| ifap | control | 0.2270722 | 0.3222 |
| ifap | all | 0.06604452  | 0.6521 |
| sCD14 | AIH | -0.5259243  | 0.0364 (*) |
| sCD14 | healthy | -0.1275773  | 0.7085 |
| sCD14 | control | -0.05840449  |  0.8015 |
| sCD14 | all | -0.2320581 | 0.1125|

**Reproduce:**

```r
source('figure6.R')

```

# citation

```
_bibtex_item_
```

