
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ConNIS <img src="./man/figures/logo.svg" alt="ConNIS" align="right" width="120"/>

<!-- badges: start -->
<!-- badges: end -->

The `ConNIS` package provides multiple methods for determining essential
genes based on insertion sites of *TraDIS* data. It also supports the
data driven determination of tuning/threshold value selection by the
implementation of a subsample based instability approach. Currently the
follwing methods are implemented:

- the *Consecutive Non-Insertion Sites* (ConNIS) method (Hanke et
  al. (2025))
- an implementation of the *Binomial* approach of the `TSAS 2.0` package
  with an additional weigh parameter
- an implementation of the *Tn5Gaps* method of the `TRANSIT` package
  with an additional weigh parameter
- the determination of essential genes by applying the geometric
  distribution

## Installation

You can install the development version of ConNIS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("moritz-hanke/ConNIS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r

library(ConNIS)
#> Loading required package: tibble
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: parallel
#> Loading required package: gmp
#> 
#> Attaching package: 'gmp'
#> The following objects are masked from 'package:base':
#> 
#>     %*%, apply, crossprod, matrix, tcrossprod

# Use the E. coli BW 25113 dataset but only the first 250 genes
truncated_ecoli <- ecoli_bw25113[1:250,]

# load the insertion sites by Goodall, 2018, but omit all insertion sites beyond
# the 250th gene
truncated_is_pos <- is_pos[is_pos <= max(truncated_ecoli$end)]

# run ConNIS with weight = 1
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 1)

results_ConNIS
#> # A tibble: 250 × 3
#>    gene  p_value weight_value
#>    <chr>   <dbl>        <dbl>
#>  1 thrL  0.446              1
#>  2 thrA  0.0827             1
#>  3 thrB  0.149              1
#>  4 thrC  0.00338            1
#>  5 yaaX  0.326              1
#>  6 yaaA  0.0102             1
#>  7 yaaJ  0.00917            1
#>  8 talB  0.0777             1
#>  9 mog   0.00893            1
#> 10 satP  0.0265             1
#> # ℹ 240 more rows
```

Using the `p.adjust()` function we use the “Bonferroni correction” for
the multiple testing problem and select the genes that are siginificant.

``` r

results_ConNIS$adjusted_p_value <-
  p.adjust(results_ConNIS$p_value, method = "bonferroni") <= 0.05

(
  results_ConNIS %>% 
    filter(adjusted_p_value) 
  )$gene
#>  [1] "dnaK" "rpsT" "ribF" "ileS" "lspA" "ispH" "dapB" "folA" "lptD" "ftsL"
#> [11] "ftsI" "murE" "murF" "mraY" "murD" "ftsW" "murG" "murC" "ddlB" "ftsQ"
#> [21] "ftsA" "ftsZ" "lpxC" "secA" "coaE" "aceF" "lpd"  "can"  "folK" "hemL"
#> [31] "erpA" "mtn"  "dapD" "map"  "rpsB" "tsf"  "pyrH" "frr"  "dxr"  "ispU"
#> [41] "cdsA" "rseP" "bamA" "lpxD" "fabZ" "lpxA" "lpxB" "rnhB" "dnaE" "accA"
#> [51] "tilS" "proS" "gloB" "ykfM" "gmhA"
```

We declare 55 genes to be essential (note, this is not the complete data
set).

Next, we re-reun ConNIS but with a smaller weight

``` r
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 0.3)

results_ConNIS$adjusted_p_value <-
  p.adjust(results_ConNIS$p_value, method = "bonferroni") <= 0.05

(
  results_ConNIS %>% 
    filter(adjusted_p_value) 
  )$gene
#>  [1] "dnaK" "rpsT" "ribF" "ileS" "lspA" "dapB" "folA" "lptD" "ftsI" "murE"
#> [11] "murF" "mraY" "murD" "ftsW" "murG" "murC" "ftsQ" "ftsA" "ftsZ" "lpxC"
#> [21] "secA" "lpd"  "can"  "erpA" "dapD" "tsf"  "pyrH" "frr"  "dxr"  "ispU"
#> [31] "cdsA" "bamA" "lpxD" "fabZ" "lpxA" "lpxB" "dnaE" "accA" "proS"
```

Only 39 genes are declared essential since smaller weights will make it
harder to label a gene by chance as ‘essential’.

Next, we compare ConNIS with Binomial, Gemeotric and Tn5Gaps all with
the same weigth.

``` r
results_Binomial <- 
  Binomial(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 0.3)

results_Binomial$adjusted_p_value <-
  p.adjust(results_Binomial$p_value, method = "bonferroni") <= 0.05

ess_genes_Binomial <- 
  (
  results_Binomial %>% 
    filter(adjusted_p_value) 
  )$gene

results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight =  0.3)

results_ConNIS$adjusted_p_value <-
  p.adjust(results_ConNIS$p_value, method = "bonferroni") <= 0.05

ess_genes_ConNIS <- 
  (
  results_ConNIS %>% 
    filter(adjusted_p_value) 
  )$gene

results_Geometric <- 
  Geometric(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight =  0.3)

results_Geometric$adjusted_p_value <-
  p.adjust(results_Geometric$p_value, method = "bonferroni") <= 0.05

ess_genes_Geometric <- 
  (
  results_Geometric %>% 
    filter(adjusted_p_value) 
  )$gene

results_Tn5Gaps <- 
  Tn5Gaps(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 0.3)

results_Tn5Gaps$adjusted_p_value <-
  p.adjust(results_Tn5Gaps$p_value, method = "bonferroni") <= 0.05

ess_genes_Tn5Gaps <- 
  (
  results_Tn5Gaps %>% 
    filter(adjusted_p_value) 
  )$gene
```

We use a Venn diagramm to compare the results:

``` r
library(ggvenn)
#> Loading required package: grid
#> Loading required package: ggplot2

ess_genes <- list(
  `Binomial` = ess_genes_Binomial,
  `ConNIS` = ess_genes_ConNIS,
  `Geometric` = ess_genes_Geometric,
  `Tn5Gaps` = ess_genes_Tn5Gaps
)

ggvenn(ess_genes)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
