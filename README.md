# AIH analysis

## notes

code depends on `metagenomeSeq` and other libraries that must be installed via either Bioconductor or CRAN. Some libraries that were used in this analysis were removed from CRAN e.g. ('biom'). You may need install these libraries (e.g. using `devtools::install_github`) manually to run these scripts.

## data

data which will be used by the scripts is located in the `data/` folder


## figure 1

![alt text](ms/version_september_17/figure1.png)

### figure 1a

alpha diversity

```
soource('figure_1a.R')
```

script generates
- `results/figure_1a_class.pdf`

which was modified manually (using pvalues from supplementary table 2) to create part a of  `ms/version_september_17/figure1.pdf`


### figure 1b

triplot

```
soource('figure_1b.R')
```

script generates
- `results/figure_1b_triplot.pdf`

which was modified to create part b of  `ms/version_september_17/figure1.pdf`

### figure 1c

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

## figure 2

![alt text](ms/version_september_17/figure2.png)

### figure 2a

PCOA

```
soource('figure_2a.R')
```

script generates
- `results/figure_2a_print.pdf`
- `results/figure_2a_label.pdf`

which was modified manually to create part a of  `ms/version_september_17/figure2.pdf`

### figure 2b

PCOA contrained by cohort

```
soource('figure_2b.R')
```

script generates
- `results/figure_2b_constrained_by_cohort_print.pdf`
- `results/figure_2b_constrained_by_cohort_label.pdf`

which was modified manually to create part b of  `ms/version_september_17/figure2.pdf`

### figure 2c

PCOA contrained by cohort

```
soource('figure_2c.R')
```

script generates
- `results/figure_2b_constrained_by_liver_print.pdf`
- `results/figure_2b_constrained_by_liver.pdf`
- `results/figure_2b_constrained_by_liver_label.pdf`

which was modified manually to create part c of  `ms/version_september_17/figure2.pdf`


## figure 3

![alt text](ms/version_september_17/figure3.png)


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

## unsorted analysis

### correlation of alpha diversity with hep. status

- [figure 5, no groups](results/figure5/figure_5_no_cat.pdf)
- [figure 5, dataset](results/figure5/data.csv)

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

### correlation of alpha diversity with markers

- [figure 6](results/figure6/figure_6.pdf)
- [figure 6 (with mean line)](results/figure6/figure_6_all.pdf)

correlation coefficients will be printed to screen (for chao1, shannon, observed)

shannon

```
> printCorrelationCoef(chao1)
control

	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "others"), ]$adiv and df.adiv[which(df.adiv$type == "others"), ]$ifap
t = -1.2159, df = 19, p-value = 0.2389
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6275839  0.1843826
sample estimates:
       cor 
-0.2686908 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "others"), ]$adiv and df.adiv[which(df.adiv$type == "others"), ]$lbp
t = 0.47462, df = 19, p-value = 0.6405
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.3392953  0.5158291
sample estimates:
      cor 
0.1082462 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "others"), ]$adiv and df.adiv[which(df.adiv$type == "others"), ]$sCD14
t = -0.96352, df = 19, p-value = 0.3474
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.5923340  0.2380273
sample estimates:
       cor 
-0.2158375 

healthy

	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "control"), ]$adiv and df.adiv[which(df.adiv$type == "control"), ]$ifap
t = -0.21226, df = 10, p-value = 0.8362
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6171537  0.5271918
sample estimates:
        cor 
-0.06697288 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "control"), ]$adiv and df.adiv[which(df.adiv$type == "control"), ]$lbp
t = -0.72758, df = 10, p-value = 0.4836
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.7071303  0.4013207
sample estimates:
       cor 
-0.2242239 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "control"), ]$adiv and df.adiv[which(df.adiv$type == "control"), ]$sCD14
t = -0.14956, df = 9, p-value = 0.8844
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6308246  0.5670197
sample estimates:
        cor 
-0.04979144 

AIH

	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "AIH"), ]$adiv and df.adiv[which(df.adiv$type == "AIH"), ]$ifap
t = 0.58046, df = 14, p-value = 0.5708
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.3705645  0.6031700
sample estimates:
      cor 
0.1533008 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "AIH"), ]$adiv and df.adiv[which(df.adiv$type == "AIH"), ]$lbp
t = 0.70243, df = 14, p-value = 0.4939
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.3425250  0.6232134
sample estimates:
      cor 
0.1845082 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "AIH"), ]$adiv and df.adiv[which(df.adiv$type == "AIH"), ]$sCD14
t = -2.4767, df = 14, p-value = 0.02664
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.82259702 -0.07744908
sample estimates:
       cor 
-0.5519635 

all

	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "AIH" | df.adiv$type == "control" |  and df.adiv[which(df.adiv$type == "AIH" | df.adiv$type == "control" |     df.adiv$type == "others"), ]$adiv and     df.adiv$type == "others"), ]$ifap
t = -2.2723, df = 47, p-value = 0.02769
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.54738273 -0.03665839
sample estimates:
       cor 
-0.3146119 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "AIH" | df.adiv$type == "control" |  and df.adiv[which(df.adiv$type == "AIH" | df.adiv$type == "control" |     df.adiv$type == "others"), ]$adiv and     df.adiv$type == "others"), ]$lbp
t = -0.54859, df = 47, p-value = 0.5859
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.3530432  0.2060530
sample estimates:
        cor 
-0.07976532 


	Pearson's product-moment correlation

data:  df.adiv[which(df.adiv$type == "AIH" | df.adiv$type == "control" |  and df.adiv[which(df.adiv$type == "AIH" | df.adiv$type == "control" |     df.adiv$type == "others"), ]$adiv and     df.adiv$type == "others"), ]$sCD14
t = -3.8275, df = 46, p-value = 0.0003892
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6805680 -0.2409882
sample estimates:
       cor 
-0.4914704 

```

