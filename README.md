
<!-- README.md is generated from README.Rmd. Please edit that file -->

# G4Micro

<!-- badges: start -->

[![R-CMD-check](https://github.com/CarlosMoraMartinez/G4Micro/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CarlosMoraMartinez/G4Micro/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

G4Micro contains helper functions used to analyze microbiome data in the
Mora-Martinez, Molina-Mendoza et al. paper.

## Installation

You can install the development version of G4Micro from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("CarlosMoraMartinez/G4Micro")
```

## Alpha Diversity Analysis

Load a phyloseq object with filtered and rarefied counts:

``` r
library(G4Micro)

data("phobj_raref")
phobj_raref
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 820 taxa and 105 samples ]
#> sample_data() Sample Data:       [ 105 samples by 119 sample variables ]
#> tax_table()   Taxonomy Table:    [ 820 taxa by 8 taxonomic ranks ]
```

A list with default options is also loaded. Copy it and modify it to use
a custom output directory.

``` r
opt <- opt_default
opt$out <- "~/test_g4micro/"
if(!dir.exists(opt$out)) dir.create(opt$out)

restoreopt <- restauraropt_mk(opt)
```

The function will allow to reset the options to the state when was
called. To do that, execute:

``` r
opt <- restoreopt(opt)
```

Calculate alpha diversity indices:

``` r
outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

alpha_indices <-  c("Observed", "Chao1", "Shannon", "InvSimpson")

divtab <- calculateAlphaDiversityTable(phseq_obj = phobj_raref, outdir = outdir, 
                                       indices = alpha_indices, name = "AlphaDiv" )

divtab %>% select(sampleID, Condition, all_of(alpha_indices)) %>% 
  head %>% kableExtra::kable()
```

| sampleID | Condition | Observed |  Chao1 |  Shannon | InvSimpson |
|:---------|:----------|---------:|-------:|---------:|-----------:|
| 1        | Control   |      398 | 398.50 | 2.657452 |   4.633405 |
| 10       | Control   |      385 | 385.00 | 2.648691 |   6.380959 |
| 100      | Control   |      329 | 329.00 | 2.982834 |   8.265564 |
| 102      | Control   |      467 | 470.75 | 2.400962 |   3.655553 |
| 103      | Control   |      447 | 447.00 | 3.509609 |  14.020145 |
| 104      | Control   |      425 | 426.00 | 3.350646 |  11.243856 |

Test statistical differences between alpha diversity indices:

``` r

alphadif <- testDiversityDifferences(divtab, alpha_indices, 
                                     groupvars = c("Condition", "Sex"), 
                                     outdir = outdir, name = "AlphaDiv_test")

alphadif %>% kableExtra::kable()
```

| variable | groups | comparison | anova_F | anova_p | t_test | wilcox_test | shapiro_normality_test | bartlett_test | levene_test | t_corrected | wilcox_corrected |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Observed | Condition | all | 17.7978164 | 0.0000530 | 0.0000534 | 0.0000116 | 0.0619656 | 0.8568824 | 0.5272713 | 0.0002357 | 0.0000521 |
| Observed | Sex | all | 0.2479607 | 0.6195762 | 0.6314333 | 0.6967562 | 0.0619656 | 0.1402254 | 0.2483335 | 0.6314333 | 0.7039865 |
| Chao1 | Condition | all | 17.5858912 | 0.0000583 | 0.0000589 | 0.0000130 | 0.0555764 | 0.8639881 | 0.5548703 | 0.0002357 | 0.0000521 |
| Chao1 | Sex | all | 0.2715406 | 0.6034188 | 0.6160029 | 0.7039865 | 0.0555764 | 0.1298843 | 0.2503562 | 0.6314333 | 0.7039865 |
| Shannon | Condition | all | 6.5589718 | 0.0118849 | 0.0081650 | 0.0094101 | 0.0316829 | 0.0450513 | 0.0448880 | 0.0163301 | 0.0188202 |
| Shannon | Sex | all | 4.6825487 | 0.0327818 | 0.0413045 | 0.0479763 | 0.0316829 | 0.0625120 | 0.0753922 | 0.0660873 | 0.0639684 |
| InvSimpson | Condition | all | 7.3268397 | 0.0079536 | 0.0063039 | 0.0002664 | 0.0000007 | 0.2403767 | 0.8755440 | 0.0163301 | 0.0007105 |
| InvSimpson | Sex | all | 2.6117387 | 0.1091344 | 0.1110623 | 0.0334241 | 0.0000007 | 0.8501868 | 0.9598410 | 0.1480830 | 0.0534785 |

Fit linear models using a main variable and several covariates:

``` r
interestvar <- "Condition"
extravars <- c("BMI", "Age", "Sex")

models <- makeLinearModelsSingleVariable(divtab, interestvar,
                                            extravars,
                                            alpha_indices,
                                            combos=1,
                                            outdir = outdir, 
                                            name = "linmodels1")
```

Now show tests for single variables:

``` r
models$single_anovas %>% select(-mod1, -mod2, nvars, Index, 
                                model, reduced_model, Df, 
                                `Pr(>F)`, padj_all) %>% 
  kableExtra::kable()
#> Warning: 'xfun::attr()' is deprecated.
#> Use 'xfun::attr2()' instead.
#> See help("Deprecated")
#> Warning: 'xfun::attr()' is deprecated.
#> Use 'xfun::attr2()' instead.
#> See help("Deprecated")
```

| nvars | Index | model | reduced_model | Df | Sum Sq | Mean Sq | F value | Pr(\>F) | padj_all | padj_bymodel |
|---:|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | Observed | Observed ~ Condition | Observed ~ 1 | 1 | 3.523265e+04 | 3.523265e+04 | 17.7978164 | 0.0000530 | 0.0004660 | 0.0001165 |
| 0 | Observed | Observed ~ BMI | Observed ~ 1 | 1 | 1.093986e+01 | 1.093986e+01 | 0.0046904 | 0.9455340 | 0.9455340 | 0.9455340 |
| 0 | Observed | Observed ~ Age | Observed ~ 1 | 1 | 3.446667e+03 | 3.446667e+03 | 1.5062743 | 0.2225062 | 0.3310526 | 0.2225062 |
| 0 | Observed | Observed ~ Sex | Observed ~ 1 | 1 | 5.743002e+02 | 5.743002e+02 | 0.2479607 | 0.6195762 | 0.7080871 | 0.6195762 |
| 0 | Chao1 | Chao1 ~ Condition | Chao1 ~ 1 | 1 | 3.511109e+04 | 3.511109e+04 | 17.5858912 | 0.0000583 | 0.0004660 | 0.0001165 |
| 0 | Chao1 | Chao1 ~ BMI | Chao1 ~ 1 | 1 | 1.156455e+01 | 1.156455e+01 | 0.0049232 | 0.9442004 | 0.9455340 | 0.9455340 |
| 0 | Chao1 | Chao1 ~ Age | Chao1 ~ 1 | 1 | 3.554813e+03 | 3.554813e+03 | 1.5436111 | 0.2169023 | 0.3310526 | 0.2225062 |
| 0 | Chao1 | Chao1 ~ Sex | Chao1 ~ 1 | 1 | 6.330390e+02 | 6.330390e+02 | 0.2715406 | 0.6034188 | 0.7080871 | 0.6195762 |
| 0 | Shannon | Shannon ~ Condition | Shannon ~ 1 | 1 | 8.804086e-01 | 8.804086e-01 | 6.5589718 | 0.0118849 | 0.0475394 | 0.0118849 |
| 0 | Shannon | Shannon ~ BMI | Shannon ~ 1 | 1 | 2.073893e-01 | 2.073893e-01 | 1.4736690 | 0.2275987 | 0.3310526 | 0.6375174 |
| 0 | Shannon | Shannon ~ Age | Shannon ~ 1 | 1 | 4.814537e-01 | 4.814537e-01 | 3.4861927 | 0.0647264 | 0.1479461 | 0.1294529 |
| 0 | Shannon | Shannon ~ Sex | Shannon ~ 1 | 1 | 6.394895e-01 | 6.394895e-01 | 4.6825487 | 0.0327818 | 0.1049018 | 0.1311273 |
| 0 | InvSimpson | InvSimpson ~ Condition | InvSimpson ~ 1 | 1 | 1.951732e+02 | 1.951732e+02 | 7.3268397 | 0.0079536 | 0.0424190 | 0.0106048 |
| 0 | InvSimpson | InvSimpson ~ BMI | InvSimpson ~ 1 | 1 | 2.833181e+01 | 2.833181e+01 | 1.0039184 | 0.3187587 | 0.4250116 | 0.6375174 |
| 0 | InvSimpson | InvSimpson ~ Age | InvSimpson ~ 1 | 1 | 9.930915e+01 | 9.930915e+01 | 3.6022252 | 0.0605020 | 0.1479461 | 0.1294529 |
| 0 | InvSimpson | InvSimpson ~ Sex | InvSimpson ~ 1 | 1 | 7.267788e+01 | 7.267788e+01 | 2.6117387 | 0.1091344 | 0.2182687 | 0.2182687 |

Show tests using covariates:

``` r
models$anovas %>% select(-mod1, -mod2, nvars, Index, 
                                model, reduced_model, Df, 
                                `Pr(>F)`, padj_all) %>% 
  kableExtra::kable()
#> Warning: 'xfun::attr()' is deprecated.
#> Use 'xfun::attr2()' instead.
#> See help("Deprecated")
#> Warning: 'xfun::attr()' is deprecated.
#> Use 'xfun::attr2()' instead.
#> See help("Deprecated")
```

| nvars | Index | model | reduced_model | Res.Df | RSS | Df | Sum of Sq | F | Pr(\>F) | padj_all | padj_bymodel |
|---:|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Observed | Observed ~ BMI + Condition | Observed ~ BMI | 101 | 235573.99217 | -1 | -3.585624e+04 | 17.953457 | 0.0000505 | 0.0002016 | 0.0001106 |
| 1 | Observed | Observed ~ Age + Condition | Observed ~ Age | 103 | 235685.29487 | -1 | -3.179639e+04 | 15.906859 | 0.0001253 | 0.0002803 | 0.0002803 |
| 1 | Observed | Observed ~ Sex + Condition | Observed ~ Sex | 103 | 238557.66170 | -1 | -3.492959e+04 | 17.496696 | 0.0000610 | 0.0002016 | 0.0001344 |
| 1 | Chao1 | Chao1 ~ BMI + Condition | Chao1 ~ BMI | 101 | 237245.79132 | -1 | -3.575986e+04 | 17.748070 | 0.0000553 | 0.0002016 | 0.0001106 |
| 1 | Chao1 | Chao1 ~ Age + Condition | Chao1 ~ Age | 103 | 237200.78334 | -1 | -3.157463e+04 | 15.662465 | 0.0001402 | 0.0002803 | 0.0002803 |
| 1 | Chao1 | Chao1 ~ Sex + Condition | Chao1 ~ Sex | 103 | 240122.55767 | -1 | -3.479060e+04 | 17.282460 | 0.0000672 | 0.0002016 | 0.0001344 |
| 1 | Shannon | Shannon ~ BMI + Condition | Shannon ~ BMI | 101 | 14.21372 | -1 | -8.394543e-01 | 6.276637 | 0.0138488 | 0.0166185 | 0.0138488 |
| 1 | Shannon | Shannon ~ Age + Condition | Shannon ~ Age | 103 | 14.22461 | -1 | -5.884919e-01 | 4.401999 | 0.0383674 | 0.0383674 | 0.0383674 |
| 1 | Shannon | Shannon ~ Sex + Condition | Shannon ~ Sex | 103 | 14.06657 | -1 | -8.227205e-01 | 6.336335 | 0.0133857 | 0.0166185 | 0.0133857 |
| 1 | InvSimpson | InvSimpson ~ BMI + Condition | InvSimpson ~ BMI | 101 | 2850.34382 | -1 | -1.943362e+02 | 7.316854 | 0.0080310 | 0.0135770 | 0.0107080 |
| 1 | InvSimpson | InvSimpson ~ Age + Condition | InvSimpson ~ Age | 103 | 2839.59009 | -1 | -1.330645e+02 | 5.014763 | 0.0273031 | 0.0297852 | 0.0364041 |
| 1 | InvSimpson | InvSimpson ~ Sex + Condition | InvSimpson ~ Sex | 103 | 2866.22137 | -1 | -1.860592e+02 | 7.080927 | 0.0090513 | 0.0135770 | 0.0120684 |

Finally, plot differences (recalculating statistical tests):

``` r
library(cowplot)

cat_vars <- c("Condition", "Sex")
num_vars <- c("BMI", "bristol_scale")

divplots <- getAlphaDiversity(phobj_raref, 
                              vars = cat_vars, 
                              qvars = num_vars,
                              opt,
                              indices = alpha_indices,
                              correct_pvalues = T, correct_pvalues_indices = T,
                              name = "alphaplots1", w = 10, h = 4)

cowplot::plot_grid(plotlist = divplots, nrow = 4)
```

<img src="man/figures/README-alpha6-1.png" width="100%" />

## Beta Diversity Analysis

PCoA on Bray-Curtis distances:

``` r
library(cowplot)

vars <- c("Condition", "Sex", "BMI")
betaplots <- makeAllPCoAs(phobj_raref, outdir,
                          method = "PCoA",
                          name = "PCoA_Bray",
                          dist_type = "bray",
                          dist_name = "Bray-Curtis",
                          vars2plot = vars,
                          extradims = 2:3,
                          labelsamples = "sampleID",
                          create_pdfs = T)

cowplot::plot_grid(plotlist = betaplots, nrow = 3)
```

<img src="man/figures/README-beta1-1.png" width="100%" />

NMDS on Bray-Curtis distances:

``` r
library(cowplot)

vars <- c("Condition", "Sex")
betaplots <- makeAllPCoAs(phobj_raref, outdir,
                          method = "NMDS",
                          name = "NMDS_Bray",
                          dist_type = "bray",
                          dist_name = "Bray-Curtis",
                          vars2plot = vars,
                          extradims = 2,
                          labelsamples = "sampleID",
                          create_pdfs = T)
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.2424031 
#> Run 1 stress 0.2412871 
#> ... New best solution
#> ... Procrustes: rmse 0.02874238  max resid 0.2059943 
#> Run 2 stress 0.2449775 
#> Run 3 stress 0.2410742 
#> ... New best solution
#> ... Procrustes: rmse 0.007181517  max resid 0.05755251 
#> Run 4 stress 0.2413643 
#> ... Procrustes: rmse 0.01549366  max resid 0.1328495 
#> Run 5 stress 0.2460339 
#> Run 6 stress 0.2424788 
#> Run 7 stress 0.2418179 
#> Run 8 stress 0.241826 
#> Run 9 stress 0.2410743 
#> ... Procrustes: rmse 0.0002787262  max resid 0.001949547 
#> ... Similar to previous best
#> Run 10 stress 0.242898 
#> Run 11 stress 0.2415093 
#> ... Procrustes: rmse 0.009435753  max resid 0.07944365 
#> Run 12 stress 0.2414201 
#> ... Procrustes: rmse 0.0129398  max resid 0.1257415 
#> Run 13 stress 0.2468357 
#> Run 14 stress 0.2467231 
#> Run 15 stress 0.242693 
#> Run 16 stress 0.2423018 
#> Run 17 stress 0.2410744 
#> ... Procrustes: rmse 0.0002642612  max resid 0.001990086 
#> ... Similar to previous best
#> Run 18 stress 0.2448093 
#> Run 19 stress 0.2447419 
#> Run 20 stress 0.2425721 
#> *** Best solution repeated 2 times

cowplot::plot_grid(plotlist = betaplots, nrow = 1)
```

<img src="man/figures/README-beta2-1.png" width="100%" />

Perform PERMANOVA with the function:

``` r
exclude_vars <- names(sample_data(phobj_raref))
exclude_vars <- exclude_vars[!exclude_vars %in% c("Condition", "Sex", "smoking_status")]
result <- makePermanova(phobj_raref,
              dist_method = "bray", 
              seed = 123, 
              exclude_vars = exclude_vars, 
              outname = "permatest") 

result %>% kableExtra::kable()
```

| variable | DF_var | DF_Residual | DF_Total | SumOfSQs_var | SumOfSQs_Residual | SumOfSQs_Total | R2_var | R2_Residual | R2_Total | F_statistic | P | perm_disp_P | perm_disp_F | perm_disp_npermuts | capscaleanova_P | capscaleanopva_F | padj | perm_disp_Padj | capscaleanova_Padj |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Condition | 1 | 103 | 104 | 0.8948514 | 16.24310 | 17.13795 | 0.0522146 | 0.9477854 | 1 | 5.6743915 | 0.0009990 | 0.0009990 | 16.924618 | 1000 | 0.0009990 | 5.2391028 | 0.0029970 | 0.0029970 | 0.0029970 |
| Sex | 1 | 103 | 104 | 0.3830715 | 16.75488 | 17.13795 | 0.0223522 | 0.9776478 | 1 | 2.3549177 | 0.0059940 | 0.0719281 | 3.406087 | 1000 | 0.0049950 | 2.2296666 | 0.0089910 | 0.1078921 | 0.0074925 |
| smoking_status | 1 | 92 | 93 | 0.1354557 | 14.70321 | 14.83867 | 0.0091286 | 0.9908714 | 1 | 0.8475647 | 0.6143856 | 0.1788212 | 2.038981 | 1000 | 0.5684316 | 0.8872324 | 0.6143856 | 0.1788212 | 0.5684316 |

Perform PERMANOVA with the function, using covariates:

``` r
permaformulas <- c(
  "braydist ~ Condition + Sex",
  "braydist ~ Condition + BMI",
  "braydist ~ Condition + Age",
  "braydist ~ Condition + Sex + BMI",
  "braydist ~ Condition + BMI + Sex + Age"
)

result <- makePermanovaFormulas(phobj_raref,
                  permaformulas,
                  dist_method = "bray", 
                  seed = 123, 
                  outname = "adonis2formulas")

result$res %>% select(model, variable, DF_var, DF_Residual, DF_Total, R2_var, R2_Residual, R2_Total, F_statistic, P, padj) %>% kableExtra::kable()
```

| model | variable | DF_var | DF_Residual | DF_Total | R2_var | R2_Residual | R2_Total | F_statistic | P | padj |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| braydist ~ Condition + Sex | Condition | 1 | 1 | 1 | 0.0522146 | 0.0522146 | 0.0522146 | 5.744691 | 0.001 | 0.0026000 |
| braydist ~ Condition + Sex | Sex | 1 | 1 | 1 | 0.0206875 | 0.0206875 | 0.0206875 | 2.276055 | 0.006 | 0.0130000 |
| braydist ~ Condition + Sex | Residual | 102 | 102 | 102 | 0.9270979 | 0.9270979 | 0.9270979 | NA | NA | NA |
| braydist ~ Condition + Sex | Total | 104 | 104 | 104 | 1.0000000 | 1.0000000 | 1.0000000 | NA | NA | NA |
| braydist ~ Condition + BMI | Condition | 1 | 1 | 1 | 0.0538639 | 0.0538639 | 0.0538639 | 5.787504 | 0.001 | 0.0026000 |
| braydist ~ Condition + BMI | BMI | 1 | 1 | 1 | 0.0154440 | 0.0154440 | 0.0154440 | 1.659407 | 0.054 | 0.0638182 |
| braydist ~ Condition + BMI | Residual | 100 | 100 | 100 | 0.9306922 | 0.9306922 | 0.9306922 | NA | NA | NA |
| braydist ~ Condition + BMI | Total | 102 | 102 | 102 | 1.0000000 | 1.0000000 | 1.0000000 | NA | NA | NA |
| braydist ~ Condition + Age | Condition | 1 | 1 | 1 | 0.0522146 | 0.0522146 | 0.0522146 | 5.694023 | 0.001 | 0.0026000 |
| braydist ~ Condition + Age | Age | 1 | 1 | 1 | 0.0124377 | 0.0124377 | 0.0124377 | 1.356340 | 0.146 | 0.1581667 |
| braydist ~ Condition + Age | Residual | 102 | 102 | 102 | 0.9353476 | 0.9353476 | 0.9353476 | NA | NA | NA |
| braydist ~ Condition + Age | Total | 104 | 104 | 104 | 1.0000000 | 1.0000000 | 1.0000000 | NA | NA | NA |
| braydist ~ Condition + Sex + BMI | Condition | 1 | 1 | 1 | 0.0538639 | 0.0538639 | 0.0538639 | 5.855260 | 0.001 | 0.0026000 |
| braydist ~ Condition + Sex + BMI | Sex | 1 | 1 | 1 | 0.0196885 | 0.0196885 | 0.0196885 | 2.140238 | 0.013 | 0.0211250 |
| braydist ~ Condition + Sex + BMI | BMI | 1 | 1 | 1 | 0.0157244 | 0.0157244 | 0.0157244 | 1.709314 | 0.049 | 0.0638182 |
| braydist ~ Condition + Sex + BMI | Residual | 99 | 99 | 99 | 0.9107233 | 0.9107233 | 0.9107233 | NA | NA | NA |
| braydist ~ Condition + Sex + BMI | Total | 102 | 102 | 102 | 1.0000000 | 1.0000000 | 1.0000000 | NA | NA | NA |
| braydist ~ Condition + BMI + Sex + Age | Condition | 1 | 1 | 1 | 0.0538639 | 0.0538639 | 0.0538639 | 5.866420 | 0.001 | 0.0026000 |
| braydist ~ Condition + BMI + Sex + Age | BMI | 1 | 1 | 1 | 0.0154440 | 0.0154440 | 0.0154440 | 1.682034 | 0.052 | 0.0638182 |
| braydist ~ Condition + BMI + Sex + Age | Sex | 1 | 1 | 1 | 0.0199689 | 0.0199689 | 0.0199689 | 2.174855 | 0.013 | 0.0211250 |
| braydist ~ Condition + BMI + Sex + Age | Age | 1 | 1 | 1 | 0.0109142 | 0.0109142 | 0.0109142 | 1.188693 | 0.274 | 0.2740000 |
| braydist ~ Condition + BMI + Sex + Age | Residual | 98 | 98 | 98 | 0.8998090 | 0.8998090 | 0.8998090 | NA | NA | NA |
| braydist ~ Condition + BMI + Sex + Age | Total | 102 | 102 | 102 | 1.0000000 | 1.0000000 | 1.0000000 | NA | NA | NA |

Model objects can also be accessed:

``` r
result$modelos[[4]]
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = as.formula(form), data = sampledf, by = adonisby, na.action = na.exclude)
#>            Df SumOfSqs      R2      F Pr(>F)    
#> Condition   1   0.9096 0.05386 5.8553  0.001 ***
#> Sex         1   0.3325 0.01969 2.1402  0.013 *  
#> BMI         1   0.2655 0.01572 1.7093  0.049 *  
#> Residual   99  15.3790 0.91072                  
#> Total     102  16.8866 1.00000                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Differential Abundance Analysis with DESeq2

First, use Sex and Age as covariates. We will use the object with the
raw counts, since normalizes data internally.

Several tables and plots will be saved to the directory.

``` r
data("phobj_filtonly")
test_vars <- c("Condition", "Sex", "Age")

result <- deseq_full_pipeline(phobj_filtonly, name = "CondSexAge", vars2deseq = test_vars, opt = opt)
#> Minfreq:  0.05 , setting minsampleswithcount to  5.25
#> converting counts to integer mode
#>   the design formula contains one or more numeric variables with integer values,
#>   specifying a model with increasing fold change for higher values.
#>   did you mean for this to be a factor? if so, first convert
#>   this variable to a factor using the factor() function
#>   the design formula contains one or more numeric variables that have mean or
#>   standard deviation larger than 5 (an arbitrary threshold to trigger this message).
#>   Including numeric variables with large mean can induce collinearity with the intercept.
#>   Users should center and scale numeric variables in the design to improve GLM convergence.
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041
#> NUM COLS:  2
#> NUM COLS:  2
#> NUM COLS:  2
#> NUM COLS:  2
#> NUM COLS:  2
#> NUM COLS:  2
#> All contrasts TRUE, intersecting Taxon listNUM COLS:  2
#> NUM COLS:  2
#> NUM COLS:  2
#> NUM COLS:  2
#> Error plotDispEsts
```

Several tables and plots will be saved to the directory.

See results table for the contrast of Condition:

``` r

taxa2plot <- result$all_contrasts$Condition_Depression_vs_Control$resdf %>% filter(!is.na(padj ) & padj < 0.01) %>% 
     arrange(padj) %>% pull(taxon) 
result$all_contrasts$Condition_Depression_vs_Control$res %>% getGTTableFromRes(taxa2plot, "DAA taxa in Depressed vs Controls")
```

<div id="joienygyjl" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#joienygyjl table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#joienygyjl thead, #joienygyjl tbody, #joienygyjl tfoot, #joienygyjl tr, #joienygyjl td, #joienygyjl th {
  border-style: none;
}
&#10;#joienygyjl p {
  margin: 0;
  padding: 0;
}
&#10;#joienygyjl .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 12px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#joienygyjl .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#joienygyjl .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#joienygyjl .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#joienygyjl .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#joienygyjl .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#joienygyjl .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#joienygyjl .gt_col_heading {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#joienygyjl .gt_column_spanner_outer {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#joienygyjl .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#joienygyjl .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#joienygyjl .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#joienygyjl .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#joienygyjl .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#joienygyjl .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#joienygyjl .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#joienygyjl .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#joienygyjl .gt_row {
  padding-top: 16px;
  padding-bottom: 16px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#joienygyjl .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#joienygyjl .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#joienygyjl .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#joienygyjl .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#joienygyjl .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#joienygyjl .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#joienygyjl .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#joienygyjl .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#joienygyjl .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#joienygyjl .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#joienygyjl .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#joienygyjl .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#joienygyjl .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#joienygyjl .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#joienygyjl .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#joienygyjl .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#joienygyjl .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#joienygyjl .gt_left {
  text-align: left;
}
&#10;#joienygyjl .gt_center {
  text-align: center;
}
&#10;#joienygyjl .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#joienygyjl .gt_font_normal {
  font-weight: normal;
}
&#10;#joienygyjl .gt_font_bold {
  font-weight: bold;
}
&#10;#joienygyjl .gt_font_italic {
  font-style: italic;
}
&#10;#joienygyjl .gt_super {
  font-size: 65%;
}
&#10;#joienygyjl .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#joienygyjl .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#joienygyjl .gt_indent_1 {
  text-indent: 5px;
}
&#10;#joienygyjl .gt_indent_2 {
  text-indent: 10px;
}
&#10;#joienygyjl .gt_indent_3 {
  text-indent: 15px;
}
&#10;#joienygyjl .gt_indent_4 {
  text-indent: 20px;
}
&#10;#joienygyjl .gt_indent_5 {
  text-indent: 25px;
}
&#10;#joienygyjl .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#joienygyjl div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>DAA taxa in Depressed vs Controls</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Taxa">Taxa</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="baseMean">baseMean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="log2FoldChange">log2FoldChange</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lfcSE">lfcSE</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="stat">stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pvalue">pvalue</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="padj">padj</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Porphyromonas somerae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1.054823e+03</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BE4B5D; color: #FFFFFF;">1.7439251</td>
<td headers="lfcSE" class="gt_row gt_right">0.3872255</td>
<td headers="stat" class="gt_row gt_right">4.503643</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2828; color: #FFFFFF;">6.679851e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CC2B2B; color: #FFFFFF;">2.605142e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Alistipes onderdonkii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">2.168480e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C24451; color: #FFFFFF;">2.2218931</td>
<td headers="lfcSE" class="gt_row gt_right">0.4206072</td>
<td headers="stat" class="gt_row gt_right">5.282584</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.273743e-07</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">6.715205e-06</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bacteroides thetaiotaomicron</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1.597751e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BB5066; color: #FFFFFF;">1.4071715</td>
<td headers="lfcSE" class="gt_row gt_right">0.3412204</td>
<td headers="stat" class="gt_row gt_right">4.123938</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CB3032; color: #FFFFFF;">3.724497e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C8383D; color: #FFFFFF;">1.162043e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bacteroides fragilis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1.454916e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #B95169; color: #FFFFFF;">1.2945258</td>
<td headers="lfcSE" class="gt_row gt_right">0.3205866</td>
<td headers="stat" class="gt_row gt_right">4.037990</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C93437; color: #FFFFFF;">5.391107e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C73A40; color: #FFFFFF;">1.298176e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bacteroides caccae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">3.212071e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C14654; color: #FFFFFF;">2.0895223</td>
<td headers="lfcSE" class="gt_row gt_right">0.3957380</td>
<td headers="stat" class="gt_row gt_right">5.280064</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.291386e-07</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">6.715205e-06</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>[Clostridium] innocuum</em></span></td>
<td headers="baseMean" class="gt_row gt_right">5.453709e+04</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #B55673; color: #FFFFFF;">0.9012218</td>
<td headers="lfcSE" class="gt_row gt_right">0.2219522</td>
<td headers="stat" class="gt_row gt_right">4.060432</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CA3335; color: #FFFFFF;">4.898192e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C73A40; color: #FFFFFF;">1.298176e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Heliomicrobium modesticaldum</em></span></td>
<td headers="baseMean" class="gt_row gt_right">2.735609e+01</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">-3.5629154</td>
<td headers="lfcSE" class="gt_row gt_right">0.9438923</td>
<td headers="stat" class="gt_row gt_right">-3.774705</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C04857; color: #FFFFFF;">1.601968e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #BA5067; color: #FFFFFF;">3.332094e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Flintibacter sp. KGMB00164</em></span></td>
<td headers="baseMean" class="gt_row gt_right">2.265668e+04</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BD4C60; color: #FFFFFF;">1.6451721</td>
<td headers="lfcSE" class="gt_row gt_right">0.2880636</td>
<td headers="stat" class="gt_row gt_right">5.711142</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.122204e-08</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.750638e-06</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Hungatella hathewayi</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1.088755e+04</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">3.8515856</td>
<td headers="lfcSE" class="gt_row gt_right">0.4153012</td>
<td headers="stat" class="gt_row gt_right">9.274198</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.789599e-20</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">5.583549e-18</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Clostridium perfringens</em></span></td>
<td headers="baseMean" class="gt_row gt_right">3.983317e+03</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BC4D61; color: #FFFFFF;">1.5902221</td>
<td headers="lfcSE" class="gt_row gt_right">0.3938914</td>
<td headers="stat" class="gt_row gt_right">4.037210</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C93437; color: #FFFFFF;">5.409066e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C73A40; color: #FFFFFF;">1.298176e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Eubacterium sp. c-25</em></span></td>
<td headers="baseMean" class="gt_row gt_right">7.560067e+03</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BE4B5D; color: #FFFFFF;">1.7303507</td>
<td headers="lfcSE" class="gt_row gt_right">0.3129395</td>
<td headers="stat" class="gt_row gt_right">5.529345</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">3.214289e-08</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">3.342861e-06</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Sellimonas intestinalis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1.854280e+04</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BC4D62; color: #FFFFFF;">1.5709083</td>
<td headers="lfcSE" class="gt_row gt_right">0.4263707</td>
<td headers="stat" class="gt_row gt_right">3.684372</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #B8536C; color: #FFFFFF;">2.292668e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #B15A7B; color: #FFFFFF;">4.378360e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Anaerotignum propionicum</em></span></td>
<td headers="baseMean" class="gt_row gt_right">8.136742e+02</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #B95169; color: #FFFFFF;">1.2984304</td>
<td headers="lfcSE" class="gt_row gt_right">0.2796630</td>
<td headers="stat" class="gt_row gt_right">4.642839</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2727; color: #FFFFFF;">3.436539e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CC2929; color: #FFFFFF;">1.531715e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Anaerostipes caccae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">2.159415e+03</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BB4E64; color: #FFFFFF;">1.4914601</td>
<td headers="lfcSE" class="gt_row gt_right">0.4271490</td>
<td headers="stat" class="gt_row gt_right">3.491663</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #8674BA; color: #FFFFFF;">4.800242e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #7D77C3; color: #FFFFFF;">7.882502e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Coprococcus eutactus</em></span></td>
<td headers="baseMean" class="gt_row gt_right">9.711471e+04</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #906FAF; color: #FFFFFF;">-1.3251483</td>
<td headers="lfcSE" class="gt_row gt_right">0.3879719</td>
<td headers="stat" class="gt_row gt_right">-3.415578</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">6.364671e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">9.928886e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Roseburia sp. NSJ-69</em></span></td>
<td headers="baseMean" class="gt_row gt_right">5.181225e+04</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #966CA8; color: #FFFFFF;">-1.0621359</td>
<td headers="lfcSE" class="gt_row gt_right">0.2694231</td>
<td headers="stat" class="gt_row gt_right">-3.942261</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C7393F; color: #FFFFFF;">8.071720e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C44049; color: #FFFFFF;">1.798840e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Roseburia intestinalis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">9.953295e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #956DA9; color: #FFFFFF;">-1.1025502</td>
<td headers="lfcSE" class="gt_row gt_right">0.3076788</td>
<td headers="stat" class="gt_row gt_right">-3.583445</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #A8628E; color: #FFFFFF;">3.390916e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #A06799; color: #FFFFFF;">5.877588e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Ruthenibacterium lactatiformans</em></span></td>
<td headers="baseMean" class="gt_row gt_right">2.349524e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #B75470; color: #FFFFFF;">1.0280344</td>
<td headers="lfcSE" class="gt_row gt_right">0.2797956</td>
<td headers="stat" class="gt_row gt_right">3.674233</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #B7546F; color: #FFFFFF;">2.385645e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #B15A7B; color: #FFFFFF;">4.378360e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Flavonifractor plautii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1.781733e+05</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #B65571; color: #FFFFFF;">0.9837113</td>
<td headers="lfcSE" class="gt_row gt_right">0.1811525</td>
<td headers="stat" class="gt_row gt_right">5.430294</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">5.626140e-08</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">4.388390e-06</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Faecalibacterium prausnitzii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">6.187147e+06</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #986BA5; color: #FFFFFF;">-0.9536573</td>
<td headers="lfcSE" class="gt_row gt_right">0.2274175</td>
<td headers="stat" class="gt_row gt_right">-4.193422</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CB2D2F; color: #FFFFFF;">2.747780e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C93539; color: #FFFFFF;">9.525637e-04</td></tr>
  </tbody>
  &#10;  
</table>
</div>

See results table for the contrast of Sex:

``` r

taxa2plot <- result$all_contrasts$Sex_Female_vs_Male$resdf %>% filter(!is.na(padj ) & padj < 0.01) %>% 
     arrange(padj) %>% pull(taxon) 
result$all_contrasts$Sex_Female_vs_Male$res %>% getGTTableFromRes(taxa2plot, "DAA taxa in Women vs Men")
```

<div id="ypsikyuwdz" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ypsikyuwdz table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#ypsikyuwdz thead, #ypsikyuwdz tbody, #ypsikyuwdz tfoot, #ypsikyuwdz tr, #ypsikyuwdz td, #ypsikyuwdz th {
  border-style: none;
}
&#10;#ypsikyuwdz p {
  margin: 0;
  padding: 0;
}
&#10;#ypsikyuwdz .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 12px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#ypsikyuwdz .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#ypsikyuwdz .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#ypsikyuwdz .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_col_heading {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#ypsikyuwdz .gt_column_spanner_outer {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#ypsikyuwdz .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#ypsikyuwdz .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#ypsikyuwdz .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#ypsikyuwdz .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#ypsikyuwdz .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#ypsikyuwdz .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#ypsikyuwdz .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#ypsikyuwdz .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#ypsikyuwdz .gt_row {
  padding-top: 16px;
  padding-bottom: 16px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#ypsikyuwdz .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ypsikyuwdz .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#ypsikyuwdz .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#ypsikyuwdz .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#ypsikyuwdz .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ypsikyuwdz .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#ypsikyuwdz .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ypsikyuwdz .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#ypsikyuwdz .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ypsikyuwdz .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#ypsikyuwdz .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ypsikyuwdz .gt_left {
  text-align: left;
}
&#10;#ypsikyuwdz .gt_center {
  text-align: center;
}
&#10;#ypsikyuwdz .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#ypsikyuwdz .gt_font_normal {
  font-weight: normal;
}
&#10;#ypsikyuwdz .gt_font_bold {
  font-weight: bold;
}
&#10;#ypsikyuwdz .gt_font_italic {
  font-style: italic;
}
&#10;#ypsikyuwdz .gt_super {
  font-size: 65%;
}
&#10;#ypsikyuwdz .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#ypsikyuwdz .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#ypsikyuwdz .gt_indent_1 {
  text-indent: 5px;
}
&#10;#ypsikyuwdz .gt_indent_2 {
  text-indent: 10px;
}
&#10;#ypsikyuwdz .gt_indent_3 {
  text-indent: 15px;
}
&#10;#ypsikyuwdz .gt_indent_4 {
  text-indent: 20px;
}
&#10;#ypsikyuwdz .gt_indent_5 {
  text-indent: 25px;
}
&#10;#ypsikyuwdz .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#ypsikyuwdz div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>DAA taxa in Women vs Men</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Taxa">Taxa</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="baseMean">baseMean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="log2FoldChange">log2FoldChange</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lfcSE">lfcSE</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="stat">stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pvalue">pvalue</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="padj">padj</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Campylobacter coli</em></span></td>
<td headers="baseMean" class="gt_row gt_right">12466.1483</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C3424D; color: #FFFFFF;">1.1638314</td>
<td headers="lfcSE" class="gt_row gt_right">0.3224306</td>
<td headers="stat" class="gt_row gt_right">3.609556</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #A56492; color: #FFFFFF;">3.067210e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #916FAE; color: #FFFFFF;">5.904379e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Campylobacter jejuni</em></span></td>
<td headers="baseMean" class="gt_row gt_right">7644.1763</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C3424D; color: #FFFFFF;">1.1797071</td>
<td headers="lfcSE" class="gt_row gt_right">0.3198959</td>
<td headers="stat" class="gt_row gt_right">3.687784</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #B45775; color: #FFFFFF;">2.262152e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #A56493; color: #FFFFFF;">4.750520e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Phocaeicola coprophilus</em></span></td>
<td headers="baseMean" class="gt_row gt_right">35206.0604</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #8574BB; color: #FFFFFF;">-1.6707527</td>
<td headers="lfcSE" class="gt_row gt_right">0.3866658</td>
<td headers="stat" class="gt_row gt_right">-4.320921</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CC2B2C; color: #FFFFFF;">1.553790e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CA3234; color: #FFFFFF;">5.982092e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Eggerthella guodeyinii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">390.4156</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C24451; color: #FFFFFF;">1.0717344</td>
<td headers="lfcSE" class="gt_row gt_right">0.2753820</td>
<td headers="stat" class="gt_row gt_right">3.891810</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C44049; color: #FFFFFF;">9.949903e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #BD4D60; color: #FFFFFF;">2.553808e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Parolsenella catena</em></span></td>
<td headers="baseMean" class="gt_row gt_right">4153.6142</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">-2.9166653</td>
<td headers="lfcSE" class="gt_row gt_right">0.4685898</td>
<td headers="stat" class="gt_row gt_right">-6.224347</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">4.835651e-10</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.117035e-07</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Olsenella timonensis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">3649.7307</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">2.1889881</td>
<td headers="lfcSE" class="gt_row gt_right">0.4900689</td>
<td headers="stat" class="gt_row gt_right">4.466695</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CC2929; color: #FFFFFF;">7.943742e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CA3133; color: #FFFFFF;">5.475345e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Actinomyces oris</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1350.1991</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C53E47; color: #FFFFFF;">1.3263722</td>
<td headers="lfcSE" class="gt_row gt_right">0.3746626</td>
<td headers="stat" class="gt_row gt_right">3.540178</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #8C71B4; color: #FFFFFF;">3.998573e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #727ACB; color: #FFFFFF;">7.105157e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Actinomyces sp. oral taxon 414</em></span></td>
<td headers="baseMean" class="gt_row gt_right">294.7989</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C3434F; color: #FFFFFF;">1.1067955</td>
<td headers="lfcSE" class="gt_row gt_right">0.3204506</td>
<td headers="stat" class="gt_row gt_right">3.453872</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">5.525991e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">8.510025e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bifidobacterium catenulatum</em></span></td>
<td headers="baseMean" class="gt_row gt_right">19169.2452</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #CC292A; color: #FFFFFF;">2.0966284</td>
<td headers="lfcSE" class="gt_row gt_right">0.4466261</td>
<td headers="stat" class="gt_row gt_right">4.694370</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2727; color: #FFFFFF;">2.674294e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CB2C2D; color: #FFFFFF;">3.088809e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Streptococcus parasanguinis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">8643.7089</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #7E76C2; color: #FFFFFF;">-1.8480696</td>
<td headers="lfcSE" class="gt_row gt_right">0.4784443</td>
<td headers="stat" class="gt_row gt_right">-3.862664</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C3424E; color: #FFFFFF;">1.121571e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #BC4D61; color: #FFFFFF;">2.590829e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Streptococcus sanguinis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1156.4990</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #8574BB; color: #FFFFFF;">-1.6682021</td>
<td headers="lfcSE" class="gt_row gt_right">0.4199901</td>
<td headers="stat" class="gt_row gt_right">-3.972003</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C73A3F; color: #FFFFFF;">7.127067e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C14755; color: #FFFFFF;">2.057941e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Clostridium perfringens</em></span></td>
<td headers="baseMean" class="gt_row gt_right">3983.3167</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #8A72B6; color: #FFFFFF;">-1.5596791</td>
<td headers="lfcSE" class="gt_row gt_right">0.3781407</td>
<td headers="stat" class="gt_row gt_right">-4.124599</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CA3133; color: #FFFFFF;">3.713805e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C63B42; color: #FFFFFF;">1.225556e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Peptacetobacter hiranonis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">154.4258</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #CA3032; color: #FFFFFF;">1.8817457</td>
<td headers="lfcSE" class="gt_row gt_right">0.4331240</td>
<td headers="stat" class="gt_row gt_right">4.344589</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CC2A2B; color: #FFFFFF;">1.395369e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CA3234; color: #FFFFFF;">5.982092e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Flavonifractor plautii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">178173.3435</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BE4A5C; color: #FFFFFF;">0.7701828</td>
<td headers="lfcSE" class="gt_row gt_right">0.1739079</td>
<td headers="stat" class="gt_row gt_right">4.428681</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CC292A; color: #FFFFFF;">9.481117e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CA3133; color: #FFFFFF;">5.475345e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Dysosmobacter sp. Marseille-Q4140</em></span></td>
<td headers="baseMean" class="gt_row gt_right">41510.1873</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BD4C5E; color: #FFFFFF;">0.7045530</td>
<td headers="lfcSE" class="gt_row gt_right">0.2020635</td>
<td headers="stat" class="gt_row gt_right">3.486790</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #637DD6; color: #FFFFFF;">4.888548e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #4882E3; color: #FFFFFF;">8.066105e-03</td></tr>
  </tbody>
  &#10;  
</table>
</div>

See results table for the contrast of Age:

``` r

taxa2plot <- result$all_contrasts$Age$resdf %>% filter(!is.na(padj ) & padj < 0.01) %>% 
     arrange(padj) %>% pull(taxon) 
result$all_contrasts$Age$res %>% getGTTableFromRes(taxa2plot, "DAA taxa with Age")
```

<div id="rkrhxkfumk" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rkrhxkfumk table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#rkrhxkfumk thead, #rkrhxkfumk tbody, #rkrhxkfumk tfoot, #rkrhxkfumk tr, #rkrhxkfumk td, #rkrhxkfumk th {
  border-style: none;
}
&#10;#rkrhxkfumk p {
  margin: 0;
  padding: 0;
}
&#10;#rkrhxkfumk .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 12px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#rkrhxkfumk .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#rkrhxkfumk .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#rkrhxkfumk .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_col_heading {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#rkrhxkfumk .gt_column_spanner_outer {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#rkrhxkfumk .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#rkrhxkfumk .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#rkrhxkfumk .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#rkrhxkfumk .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#rkrhxkfumk .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#rkrhxkfumk .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#rkrhxkfumk .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#rkrhxkfumk .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#rkrhxkfumk .gt_row {
  padding-top: 16px;
  padding-bottom: 16px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#rkrhxkfumk .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rkrhxkfumk .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#rkrhxkfumk .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#rkrhxkfumk .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#rkrhxkfumk .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rkrhxkfumk .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#rkrhxkfumk .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rkrhxkfumk .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#rkrhxkfumk .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rkrhxkfumk .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#rkrhxkfumk .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rkrhxkfumk .gt_left {
  text-align: left;
}
&#10;#rkrhxkfumk .gt_center {
  text-align: center;
}
&#10;#rkrhxkfumk .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#rkrhxkfumk .gt_font_normal {
  font-weight: normal;
}
&#10;#rkrhxkfumk .gt_font_bold {
  font-weight: bold;
}
&#10;#rkrhxkfumk .gt_font_italic {
  font-style: italic;
}
&#10;#rkrhxkfumk .gt_super {
  font-size: 65%;
}
&#10;#rkrhxkfumk .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#rkrhxkfumk .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#rkrhxkfumk .gt_indent_1 {
  text-indent: 5px;
}
&#10;#rkrhxkfumk .gt_indent_2 {
  text-indent: 10px;
}
&#10;#rkrhxkfumk .gt_indent_3 {
  text-indent: 15px;
}
&#10;#rkrhxkfumk .gt_indent_4 {
  text-indent: 20px;
}
&#10;#rkrhxkfumk .gt_indent_5 {
  text-indent: 25px;
}
&#10;#rkrhxkfumk .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#rkrhxkfumk div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>DAA taxa with Age</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Taxa">Taxa</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="baseMean">baseMean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="log2FoldChange">log2FoldChange</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lfcSE">lfcSE</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="stat">stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pvalue">pvalue</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="padj">padj</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Enterobacter cloacae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">55.32436</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">0.78327702</td>
<td headers="lfcSE" class="gt_row gt_right">0.11267862</td>
<td headers="stat" class="gt_row gt_right">6.951426</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">3.616133e-12</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.471766e-09</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Parabacteroides goldsteinii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">20872.76272</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #7579C9; color: #FFFFFF;">0.07264734</td>
<td headers="lfcSE" class="gt_row gt_right">0.01667041</td>
<td headers="stat" class="gt_row gt_right">4.357861</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C8373B; color: #FFFFFF;">1.313399e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C24450; color: #FFFFFF;">1.375706e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Prevotella buccalis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">3763.69919</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #1E86EE; color: #FFFFFF;">-0.07566150</td>
<td headers="lfcSE" class="gt_row gt_right">0.01779742</td>
<td headers="stat" class="gt_row gt_right">-4.251262</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C53F47; color: #FFFFFF;">2.125689e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C24450; color: #FFFFFF;">1.375706e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bacteroides fragilis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">145491.59199</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #3A84E8; color: #FFFFFF;">-0.05079489</td>
<td headers="lfcSE" class="gt_row gt_right">0.01201618</td>
<td headers="stat" class="gt_row gt_right">-4.227209</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C4414B; color: #FFFFFF;">2.366079e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C24450; color: #FFFFFF;">1.375706e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Schaalia turicensis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">213.32202</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #3385EA; color: #FFFFFF;">-0.05896909</td>
<td headers="lfcSE" class="gt_row gt_right">0.01536695</td>
<td headers="stat" class="gt_row gt_right">-3.837398</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">1.243448e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">6.326043e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Streptococcus canis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">689.89313</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #707ACD; color: #FFFFFF;">0.05830204</td>
<td headers="lfcSE" class="gt_row gt_right">0.01377455</td>
<td headers="stat" class="gt_row gt_right">4.232590</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C4404A; color: #FFFFFF;">2.310159e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C24450; color: #FFFFFF;">1.375706e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Clostridium taeniosporum</em></span></td>
<td headers="baseMean" class="gt_row gt_right">12.62439</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #9D699F; color: #FFFFFF;">0.25130604</td>
<td headers="lfcSE" class="gt_row gt_right">0.05932096</td>
<td headers="stat" class="gt_row gt_right">4.236379</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C4404A; color: #FFFFFF;">2.271537e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C24450; color: #FFFFFF;">1.375706e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Sellimonas intestinalis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">18542.80200</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">-0.07650053</td>
<td headers="lfcSE" class="gt_row gt_right">0.01598123</td>
<td headers="stat" class="gt_row gt_right">-4.786899</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CC2829; color: #FFFFFF;">1.693781e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CB2F31; color: #FFFFFF;">3.446844e-04</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Make heatmap:

``` r

taxa2plot <- result$all_contrasts$Condition_Depression_vs_Control$resdf %>% filter(!is.na(padj ) & padj < 0.01) %>% 
     arrange(padj) %>% pull(taxon) 

makeHeatmap(result$all_contrasts$Condition_Depression_vs_Control$resdf, 
            result$dds, 
            result$vst_counts_df, 
            c("Condition", "Sex", "Age"),
            opt, 
            name ="test_heatmap",
            logscale = F, 
            ptype="padj", 
            trim_values = TRUE, 
            taxalist=taxa2plot, 
            max_hm_h=10, max_hm_w=12)
#> NUM COLS:  2
```

<img src="man/figures/README-deseq5-1.png" width="100%" />

Also, different contrasts can be compared:

``` r

mainContrast <- result$all_contrasts$Condition_Depression_vs_Control
contrastlist2 <- list(
  result$all_contrasts$Sex_Female_vs_Male$resdf,
  result$all_contrasts$Age$resdf
) %>% lapply(\(x)return(list(resdf=x)))

names(contrastlist2) <- c( "Sex", "Age")

contrast_names_pretty <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))

compareLFCContrats2(contrastlist2, mainContrast, 
                   contrast_names_pretty, "Depression_vs_Control", 
                   plim_select= 0.001, plim_plot=0.05,
                   name2remove = "xxx",
                   resdfname="resdf", 
                   outdir = opt$out, 
                   name="LFC_Comparison_AgeAndBMI_allCombos_p05", 
                   w=12, h=12, scale_mode = "free")
```

<img src="man/figures/README-deseq6-1.png" width="100%" />

``` r

restoreopt <- restauraropt_mk(opt)
```

## DESeq2 analysis with interaction

Instead of using several covariates in an additive way, one may specify
more complex formulas with interactions:

``` r

phobj_filtered <- subset_samples(
  phobj_filtonly,
  !is.na(Condition) &
  !is.na(Sex) &
  !is.na(BMI) &
  !is.na(Age)
)
result_int <- getDeseqResults(phobj_filtered, opt = opt, name = "CBSA_interaction", 
                          variables=NULL, formula="~ Condition + BMI + Sex * Age")
#> converting counts to integer mode
#>   the design formula contains one or more numeric variables with integer values,
#>   specifying a model with increasing fold change for higher values.
#>   did you mean for this to be a factor? if so, first convert
#>   this variable to a factor using the factor() function
#>   the design formula contains one or more numeric variables that have mean or
#>   standard deviation larger than 5 (an arbitrary threshold to trigger this message).
#>   Including numeric variables with large mean can induce collinearity with the intercept.
#>   Users should center and scale numeric variables in the design to improve GLM convergence.
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> 7 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895
#> using 'ashr' for LFC shrinkage. If used in published research, please cite:
#>     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
#>     https://doi.org/10.1093/biostatistics/kxw041




taxa2plot <- result_int$all_contrasts$Condition_Depression_vs_Control$resdf %>% 
  filter(!is.na(padj ) & padj < 0.01) %>% 
  arrange(padj) %>% pull(taxon) 
result_int$all_contrasts$Condition_Depression_vs_Control$res %>% 
  getGTTableFromRes(taxa2plot, "DAA between Depr and C controlling for  BMI + Sex * Age")
```

<div id="obdyibudeb" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#obdyibudeb table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#obdyibudeb thead, #obdyibudeb tbody, #obdyibudeb tfoot, #obdyibudeb tr, #obdyibudeb td, #obdyibudeb th {
  border-style: none;
}
&#10;#obdyibudeb p {
  margin: 0;
  padding: 0;
}
&#10;#obdyibudeb .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 12px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#obdyibudeb .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#obdyibudeb .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#obdyibudeb .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_col_heading {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#obdyibudeb .gt_column_spanner_outer {
  color: #333333;
  background-color: #D3D3D3;
  font-size: 16px;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#obdyibudeb .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#obdyibudeb .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#obdyibudeb .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#obdyibudeb .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#obdyibudeb .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#obdyibudeb .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#obdyibudeb .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#obdyibudeb .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#obdyibudeb .gt_row {
  padding-top: 16px;
  padding-bottom: 16px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#obdyibudeb .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#obdyibudeb .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#obdyibudeb .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#obdyibudeb .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#obdyibudeb .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#obdyibudeb .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#obdyibudeb .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#obdyibudeb .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#obdyibudeb .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#obdyibudeb .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#obdyibudeb .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#obdyibudeb .gt_left {
  text-align: left;
}
&#10;#obdyibudeb .gt_center {
  text-align: center;
}
&#10;#obdyibudeb .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#obdyibudeb .gt_font_normal {
  font-weight: normal;
}
&#10;#obdyibudeb .gt_font_bold {
  font-weight: bold;
}
&#10;#obdyibudeb .gt_font_italic {
  font-style: italic;
}
&#10;#obdyibudeb .gt_super {
  font-size: 65%;
}
&#10;#obdyibudeb .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#obdyibudeb .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#obdyibudeb .gt_indent_1 {
  text-indent: 5px;
}
&#10;#obdyibudeb .gt_indent_2 {
  text-indent: 10px;
}
&#10;#obdyibudeb .gt_indent_3 {
  text-indent: 15px;
}
&#10;#obdyibudeb .gt_indent_4 {
  text-indent: 20px;
}
&#10;#obdyibudeb .gt_indent_5 {
  text-indent: 25px;
}
&#10;#obdyibudeb .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#obdyibudeb div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>DAA between Depr and C controlling for  BMI + Sex * Age</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Taxa">Taxa</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="baseMean">baseMean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="log2FoldChange">log2FoldChange</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lfcSE">lfcSE</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="stat">stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pvalue">pvalue</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="padj">padj</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Porphyromonas somerae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1061.7190</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C4404A; color: #FFFFFF;">1.6760097</td>
<td headers="lfcSE" class="gt_row gt_right">0.4076235</td>
<td headers="stat" class="gt_row gt_right">4.111661</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C63C44; color: #FFFFFF;">3.928234e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #BF4A5A; color: #FFFFFF;">2.051411e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Alistipes onderdonkii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">216103.2893</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">2.5447594</td>
<td headers="lfcSE" class="gt_row gt_right">0.4271729</td>
<td headers="stat" class="gt_row gt_right">5.957212</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">2.565771e-09</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">6.029563e-07</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bacteroides thetaiotaomicron</em></span></td>
<td headers="baseMean" class="gt_row gt_right">161575.7018</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C04756; color: #FFFFFF;">1.3776995</td>
<td headers="lfcSE" class="gt_row gt_right">0.3640082</td>
<td headers="stat" class="gt_row gt_right">3.784804</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #A06799; color: #FFFFFF;">1.538299e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #8D71B3; color: #FFFFFF;">5.464936e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Bacteroides caccae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">325516.1153</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C63D45; color: #FFFFFF;">1.8056324</td>
<td headers="lfcSE" class="gt_row gt_right">0.4121648</td>
<td headers="stat" class="gt_row gt_right">4.380850</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CB2E2F; color: #FFFFFF;">1.182171e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C8363B; color: #FFFFFF;">7.937432e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Actinomyces naeslundii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">304.8557</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #717ACC; color: #FFFFFF;">-1.4346297</td>
<td headers="lfcSE" class="gt_row gt_right">0.3822156</td>
<td headers="stat" class="gt_row gt_right">-3.753456</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #956DA9; color: #FFFFFF;">1.744129e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #8D71B3; color: #FFFFFF;">5.464936e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Ligilactobacillus ruminis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">9899.6459</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">-2.2004059</td>
<td headers="lfcSE" class="gt_row gt_right">0.4782790</td>
<td headers="stat" class="gt_row gt_right">-4.600674</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CC2929; color: #FFFFFF;">4.211254e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CB2D2F; color: #FFFFFF;">3.298816e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Flintibacter sp. KGMB00164</em></span></td>
<td headers="baseMean" class="gt_row gt_right">22827.2630</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C4414B; color: #FFFFFF;">1.6581181</td>
<td headers="lfcSE" class="gt_row gt_right">0.3083450</td>
<td headers="stat" class="gt_row gt_right">5.377476</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">7.553725e-08</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.183417e-05</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Clostridium perfringens</em></span></td>
<td headers="baseMean" class="gt_row gt_right">4030.6815</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C24450; color: #FFFFFF;">1.5279191</td>
<td headers="lfcSE" class="gt_row gt_right">0.4183177</td>
<td headers="stat" class="gt_row gt_right">3.652533</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">2.596667e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #1C86EE; color: #FFFFFF;">7.627709e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Eubacterium sp. c-25</em></span></td>
<td headers="baseMean" class="gt_row gt_right">7621.5185</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C53F48; color: #FFFFFF;">1.7202547</td>
<td headers="lfcSE" class="gt_row gt_right">0.3359070</td>
<td headers="stat" class="gt_row gt_right">5.121223</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">3.035603e-07</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2727; color: #FFFFFF;">3.566833e-05</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Anaerotignum propionicum</em></span></td>
<td headers="baseMean" class="gt_row gt_right">816.5338</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BF4959; color: #FFFFFF;">1.3033137</td>
<td headers="lfcSE" class="gt_row gt_right">0.3003106</td>
<td headers="stat" class="gt_row gt_right">4.339885</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CB2F31; color: #FFFFFF;">1.425573e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #C8373C; color: #FFFFFF;">8.375241e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Enterocloster bolteae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">141754.6559</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C9363A; color: #FFFFFF;">2.0737379</td>
<td headers="lfcSE" class="gt_row gt_right">0.3447668</td>
<td headers="stat" class="gt_row gt_right">6.014901</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">1.799968e-09</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CD2626; color: #FFFFFF;">6.029563e-07</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Anaerostipes caccae</em></span></td>
<td headers="baseMean" class="gt_row gt_right">2177.6466</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #C4414B; color: #FFFFFF;">1.6429036</td>
<td headers="lfcSE" class="gt_row gt_right">0.4160806</td>
<td headers="stat" class="gt_row gt_right">3.948523</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #BD4D61; color: #FFFFFF;">7.863492e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #B15A7B; color: #FFFFFF;">3.359855e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Roseburia sp. NSJ-69</em></span></td>
<td headers="baseMean" class="gt_row gt_right">52433.4304</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #8175BF; color: #FFFFFF;">-1.1334449</td>
<td headers="lfcSE" class="gt_row gt_right">0.2816556</td>
<td headers="stat" class="gt_row gt_right">-4.024223</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #C24451; color: #FFFFFF;">5.716365e-05</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #B9526A; color: #FFFFFF;">2.686692e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Roseburia intestinalis</em></span></td>
<td headers="baseMean" class="gt_row gt_right">1006970.8929</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #7D77C2; color: #FFFFFF;">-1.2174001</td>
<td headers="lfcSE" class="gt_row gt_right">0.3210338</td>
<td headers="stat" class="gt_row gt_right">-3.792125</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #A26596; color: #FFFFFF;">1.493640e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #8D71B3; color: #FFFFFF;">5.464936e-03</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Flavonifractor plautii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">179275.8932</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #BA5067; color: #FFFFFF;">0.9571934</td>
<td headers="lfcSE" class="gt_row gt_right">0.1974312</td>
<td headers="stat" class="gt_row gt_right">4.848237</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #CD2727; color: #FFFFFF;">1.245634e-06</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #CC2929; color: #FFFFFF;">1.170896e-04</td></tr>
    <tr><td headers="Taxa" class="gt_row gt_left"><span class='gt_from_md'><em>Faecalibacterium prausnitzii</em></span></td>
<td headers="baseMean" class="gt_row gt_right">6212424.2010</td>
<td headers="log2FoldChange" class="gt_row gt_right" style="background-color: #8A72B5; color: #FFFFFF;">-0.9157742</td>
<td headers="lfcSE" class="gt_row gt_right">0.2433100</td>
<td headers="stat" class="gt_row gt_right">-3.763817</td>
<td headers="pvalue" class="gt_row gt_right" style="background-color: #996BA4; color: #FFFFFF;">1.673393e-04</td>
<td headers="padj" class="gt_row gt_right" style="background-color: #8D71B3; color: #FFFFFF;">5.464936e-03</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r

restoreopt <- restauraropt_mk(opt)
```
