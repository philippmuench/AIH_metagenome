# AIH paper

Specific alterations of the intestinal microbiome in autoimmune hepatitis show partially liver disease-specific patterns being also related to the extent of liver parenchymatous tissue remodeling

# Table of Contents  
[Usuage](#usuage)  
[Figure 1 (alpha diversity)](#figure-1)  
[Figure 2 (class level visualization)](#figure-2)  
[Citation](#citation)

# Usuage

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
```

# Figures
## Figure 1
**Description:** 
alpha diversity plots for shannon, chao1 and observed OTU
others group: NASH, alkohol, others, ALD, HCV

**File:**
- [figure1 (observed OTUs index)](results/figure1/figure_1_index_observed.pdf)
- [figure1 (chao1 index)](results/figure1/figure_1_index_chao1.pdf)
- [figure1 (Shannon index)](results/figure1/figure_1_index_shannon.pdf)

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
Applied filtering: `removeRelAb()` with 0.01 (genera with rel. abundance < 1% are removed)

**File:**
- [figure2 (merged)](results/figure2/figure_2.pdf)
- [figure2 (bar plot)](results/figure2/figure_2_bar.pdf)
- [figure2 (difference)](results/figure2/figure_2_diff.pdf)
- [figure2 (triplot)](results/figure2/figure_2_tri.pdf)

**Reproduce:**

```r
source('figure2.R')
```

# citation

```
_bibtex_item_
```

