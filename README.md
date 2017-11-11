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

### healthy vs. non-AIH hep. control

| Class*                	| logFC 	| t     	| P value  	| FDR      	| B     	| sig. Level 	|
|-----------------------	|-------	|-------	|----------	|----------	|-------	|------------	|
| RF3                   	| 1.47  	| 3.45  	| 0.003716 	| 0.027254 	| -1.92 	|  *         	|
| Opitutae              	| 1.18  	| 1.96  	| 0.083799 	| 0.307263 	| -4.45 	| n.s.       	|
| Deferribacteres       	| 1.03  	| 1.55  	| 0.142350 	| 0.377765 	| -4.80 	| n.s.       	|
| [Lentisphaeria]       	| 0.72  	| 0.87  	| 0.402788 	| 0.579287 	| -5.98 	| n.s.       	|
| Erysipelotrichi       	| 0.37  	| 1.10  	| 0.277168 	| 0.435550 	| -6.44 	| n.s.       	|
| Clostridia            	| 0.33  	| 2.32  	| 0.024737 	| 0.136056 	| -4.48 	| n.s.       	|
| 4C0d-2                	| 0.15  	| 0.30  	| 0.769840 	| 0.865654 	| -6.62 	| n.s.       	|
| Mollicutes            	| 0.11  	| 0.13  	| 0.899250 	| 0.905493 	| -6.62 	| n.s.       	|
| Deltaproteobacteria   	| -0.05 	| -0.12 	| 0.905493 	| 0.905493 	| -6.85 	| n.s.       	|
| Bacilli               	| -0.09 	| -0.27 	| 0.786958 	| 0.865654 	| -7.03 	| n.s.       	|
| Bacteroidia           	| -0.11 	| -0.57 	| 0.572722 	| 0.699994 	| -6.90 	| n.s.       	|
| Alphaproteobacteria   	| -0.18 	| -0.81 	| 0.421300 	| 0.579287 	| -6.73 	| n.s.       	|
| Betaproteobacteria    	| -0.25 	| -1.10 	| 0.276335 	| 0.435550 	| -6.44 	| n.s.       	|
| Verrucomicrobiae      	| -0.36 	| -0.60 	| 0.554026 	| 0.699994 	| -6.71 	| n.s.       	|
| Coriobacteriia        	| -0.43 	| -1.29 	| 0.204429 	| 0.408857 	| -6.20 	| n.s.       	|
| Flavobacteriia        	| -0.49 	| -1.17 	| 0.257856 	| 0.435550 	| -5.86 	| n.s.       	|
| Chloroplast           	| -0.50 	| -1.42 	| 0.166344 	| 0.377765 	| -5.76 	| n.s.       	|
| Actinobacteria        	| -0.62 	| -1.45 	| 0.154523 	| 0.377765 	| -6.02 	| n.s.       	|
| Epsilonproteobacteria 	| -0.93 	| -1.42 	| 0.171711 	| 0.377765 	| -5.54 	| n.s.       	|
| Gammaproteobacteria   	| -1.38 	| -4.05 	| 0.000203 	| 0.002229 	| 0.01  	|  * *       	|
| Deinococci            	| -1.87 	| -2.07 	| 0.073636 	| 0.307263 	| -4.22 	| n.s.       	|
| Fusobacteriia         	| -2.33 	| -6.53 	| 0.000001 	| 0.000022 	| 5.67  	|  * * *     	|

### AIH vs. healthy

| Class*                	| logFC 	| t     	| P value 	| FDR    	| B     	| sig. Level 	|
|-----------------------	|-------	|-------	|---------	|--------	|-------	|------------	|
| Fusobacteriia         	| 1.67  	| 3.39  	| 0.0027  	| 0.0314 	| -1.57 	|  *         	|
| Gammaproteobacteria   	| 1.51  	| 4.28  	| 0.0001  	| 0.0027 	| 1.04  	|  * *       	|
| Synergistia           	| 0.76  	| 0.81  	| 0.4371  	| 0.7733 	| -5.16 	| n.s.       	|
| Verrucomicrobiae      	| 0.44  	| 0.83  	| 0.4101  	| 0.7733 	| -6.00 	| n.s.       	|
| Actinobacteria        	| 0.41  	| 1.08  	| 0.2876  	| 0.7158 	| -5.87 	| n.s.       	|
| Bacteroidia           	| 0.29  	| 1.39  	| 0.1728  	| 0.6025 	| -5.50 	| n.s.       	|
| Coriobacteriia        	| 0.29  	| 1.03  	| 0.3112  	| 0.7158 	| -5.91 	| n.s.       	|
| Betaproteobacteria    	| 0.24  	| 0.94  	| 0.3523  	| 0.7365 	| -5.99 	| n.s.       	|
| Erysipelotrichi       	| 0.23  	| 0.72  	| 0.4743  	| 0.7791 	| -6.18 	| n.s.       	|
| RF3                   	| 0.21  	| 0.41  	| 0.6899  	| 0.8755 	| -5.71 	| n.s.       	|
| Bacilli               	| 0.16  	| 0.56  	| 0.5804  	| 0.8343 	| -6.29 	| n.s.       	|
| Deferribacteres       	| 0.14  	| 0.28  	| 0.7846  	| 0.8755 	| -5.90 	| n.s.       	|
| 4C0d-2                	| 0.13  	| 0.26  	| 0.7980  	| 0.8755 	| -6.06 	| n.s.       	|
| Deinococci            	| 0.11  	| 0.15  	| 0.8809  	| 0.9210 	| -5.61 	| n.s.       	|
| Chloroplast           	| 0.11  	| 0.26  	| 0.7993  	| 0.8755 	| -6.17 	| n.s.       	|
| Epsilonproteobacteria 	| 0.05  	| 0.08  	| 0.9404  	| 0.9404 	| -5.98 	| n.s.       	|
| Clostridia            	| -0.25 	| -1.32 	| 0.1937  	| 0.6025 	| -5.59 	| n.s.       	|
| Flavobacteriia        	| -0.25 	| -0.57 	| 0.5737  	| 0.8343 	| -5.86 	| n.s.       	|
| Mollicutes            	| -0.28 	| -0.38 	| 0.7050  	| 0.8755 	| -6.12 	| n.s.       	|
| Alphaproteobacteria   	| -0.37 	| -1.56 	| 0.1264  	| 0.5814 	| -5.26 	| n.s.       	|
| Deltaproteobacteria   	| -0.45 	| -1.28 	| 0.2096  	| 0.6025 	| -5.55 	| n.s.       	|
| [Lentisphaeria]       	| -1.27 	| -1.83 	| 0.0861  	| 0.4953 	| -4.37 	| n.s.       	|
| Opitutae              	| -1.52 	| -1.99 	| 0.0719  	| 0.4953 	| -4.05 	| n.s.       	|

### AIH vs. non-AIH hep. control

| Class*                	| logFC 	| t     	| P value 	| FDR  	| B     	| sig. Level 	|
|-----------------------	|-------	|-------	|---------	|------	|-------	|------------	|
| RF3                   	| 1.68  	| 2.34  	| 0.04    	| 0.51 	| -4.57 	| n.s.       	|
| Deferribacteres       	| 1.17  	| 0.85  	| 0.42    	| 0.84 	| -4.60 	| n.s.       	|
| Erysipelotrichi       	| 0.60  	| 1.54  	| 0.13    	| 0.63 	| -4.54 	| n.s.       	|
| 4C0d-2                	| 0.28  	| 0.36  	| 0.72    	| 0.92 	| -4.61 	| n.s.       	|
| Bacteroidia           	| 0.18  	| 0.93  	| 0.36    	| 0.84 	| -4.60 	| n.s.       	|
| Gammaproteobacteria   	| 0.14  	| 0.37  	| 0.71    	| 0.92 	| -4.63 	| n.s.       	|
| Clostridia            	| 0.08  	| 0.48  	| 0.64    	| 0.92 	| -4.63 	| n.s.       	|
| Verrucomicrobiae      	| 0.07  	| 0.10  	| 0.92    	| 0.96 	| -4.63 	| n.s.       	|
| Bacilli               	| 0.07  	| 0.17  	| 0.86    	| 0.94 	| -4.64 	| n.s.       	|
| Betaproteobacteria    	| 0.00  	| -0.02 	| 0.98    	| 0.98 	| -4.64 	| n.s.       	|
| Coriobacteriia        	| -0.14 	| -0.44 	| 0.66    	| 0.92 	| -4.63 	| n.s.       	|
| Mollicutes            	| -0.17 	| -0.18 	| 0.86    	| 0.94 	| -4.61 	| n.s.       	|
| Actinobacteria        	| -0.21 	| -0.53 	| 0.60    	| 0.92 	| -4.63 	| n.s.       	|
| Elusimicrobia         	| -0.28 	| -0.19 	| 0.86    	| 0.94 	| -4.60 	| n.s.       	|
| Opitutae              	| -0.34 	| -0.47 	| 0.65    	| 0.92 	| -4.60 	| n.s.       	|
| Chloroplast           	| -0.39 	| -0.76 	| 0.45    	| 0.84 	| -4.61 	| n.s.       	|
| Deltaproteobacteria   	| -0.50 	| -1.07 	| 0.30    	| 0.84 	| -4.59 	| n.s.       	|
| [Lentisphaeria]       	| -0.56 	| -0.79 	| 0.44    	| 0.84 	| -4.60 	| n.s.       	|
| Alphaproteobacteria   	| -0.56 	| -2.17 	| 0.04    	| 0.51 	| -4.45 	| n.s.       	|
| TM7-3                 	| -0.58 	| -0.85 	| 0.42    	| 0.84 	| -4.60 	| n.s.       	|
| Fusobacteriia         	| -0.66 	| -1.29 	| 0.21    	| 0.72 	| -4.58 	| n.s.       	|
| Flavobacteriia        	| -0.74 	| -1.48 	| 0.16    	| 0.64 	| -4.58 	| n.s.       	|
| Epsilonproteobacteria 	| -0.88 	| -1.69 	| 0.10    	| 0.63 	| -4.55 	| n.s.       	|
| Deinococci            	| -1.76 	| -1.69 	| 0.13    	| 0.63 	| -4.59 	| n.s.       	|


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

