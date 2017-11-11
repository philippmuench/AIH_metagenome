
# AIH analysis

## notes

code depends on `metagenomeSeq` and other libraries that must be installed via either Bioconductor or CRAN. Some libraries that were used in this analysis were removed from CRAN e.g. ('biom'). You may need install these libraries (e.g. using `devtools::install_github`) manually to run these scripts.

## figure 1a

alpha diversity

```
soource('figure_1a.R')
```

script generates
- `results/figure_1a_class.pdf`

which was modified manually (using pvalues from supplementary table 2) to create part a of  `ms/version_september_17/figure1.pdf`




## figure 1c

alpha diversity

```
soource('figure_1c.R')
```

script generates
- `results/figure_1c_chao1.pdf`
- test statistics gets printed to screen:

```
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = value ~ type, data = df.adiv)

$type
                     diff       lwr        upr     p adj
control-AIH    -0.2994906 -1.252769 0.65378825 0.7308472
others-AIH     -0.8918631 -1.806463 0.02273673 0.0574292
others-control -0.5923725 -1.463349 0.27860441 0.2384276


```

which was modified manually to create part c of  `ms/version_september_17/figure1.pdf`

## figure 3

analysis on OTU level

```
soource('figure_3.R')
```

script generates
- `results/figure_3_triplot.pdf`
- `results/figure_3_aih_db.csv`
- `results/figure_3_control_db.csv`
- `results/figure_3_healthy_db.csv`

which was modified manually to create `ms/version_september_17/figure3.pdf`

## supplementary table 2

statistics on class level

```
soource('table_s2.R')
```

script generates
- `results/table_s2_class_aih_vs_control.tsv`
- `table_s2_class_aih_vs_helathy.tsv`
- `table_s2_class_healthy_vs_control.tsv`

which was modified manually to create `ms/version_september_17/Table_S2.xlsx`

## supplementary table 3

statistics on family level

```
soource('table_s3.R')
```

script generates
- `results/table_s3_family_aih_vs_control.tsv`
- `table_s3_family_aih_vs_helathy.tsv`
- `table_s3_family_healthy_vs_control.tsv`

which was modified manually to create `ms/version_september_17/Table_S3.xlsx`

## supplementary table 4

statistics on genus level

```
soource('table_s4.R')
```

script generates
- `results/table_s4_genus_aih_vs_control.tsv`
- `table_s3_genus_aih_vs_helathy.tsv`
- `table_s3_genus_healthy_vs_control.tsv`

which was modified manually to create `ms/version_september_17/Table_S4.xlsx`



## supplementary figure 1

```
soource('figure_s1.R')
```

script generates
- `results/figure_s1a_observed.pdf`
- `results/figure_s1b_shannon.pdf`

which was modified manually to create part c of  `ms/version_september_17/figure_s3.pdf`



## supplementary figure 2

```
soource('figure_s2.R')
```

script generates
- `results/figure_s2.pdf`

which was modified manually to create part c of  `ms/version_september_17/figure_s2.pdf`


## supplementary figure 3

```
soource('figure_s3.R')
```

script generates
- `results/figure_s3.pdf`

which was modified manually to create part c of  `ms/version_september_17/figure_s3.pdf`
