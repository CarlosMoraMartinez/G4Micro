
<!-- README.md is generated from README.Rmd. Please edit that file -->

# G4Micro

<!-- badges: start -->

[![R-CMD-check](https://github.com/CarlosMoraMartinez/G4Micro/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CarlosMoraMartinez/G4Micro/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

G4Micro contains functions used to analyze microbiome data in the
Mora-Martinez, Molina-Mendoza et al. paper.

## Installation

You can install the development version of G4Micro from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("CarlosMoraMartinez/G4Micro")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(G4Micro)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: forcats
#> Loading required package: ggplot2
#> Loading required package: phyloseq
#> Loading required package: purrr
#> Loading required package: readr
#> Loading required package: stringr
#> Loading required package: tibble
#> Loading required package: tidyr
#> Warning: replacing previous import 'magrittr::set_names' by 'purrr::set_names'
#> when loading 'G4Micro'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'G4Micro'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'G4Micro'
#> Warning: replacing previous import 'magrittr::extract' by 'tidyr::extract' when
#> loading 'G4Micro'
#> Warning: replacing previous import 'stats::filter' by 'dplyr::filter' when
#> loading 'G4Micro'
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
