
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

# Use the E. coli BW 25113 dataset but only the first 100 genes
truncated_ecoli <- ecoli_bw25113[70:100,]

# load the insertion sites by Goodall, 2018, but omit all insertion sites that are
# within the truncated_ecoli
truncated_is_pos <- is_pos[is_pos <= max(truncated_ecoli$end) &
                           is_pos >= min(truncated_ecoli$start)]

# run ConNIS with weight = 1
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 1)

results_ConNIS
#> # A tibble: 31 × 3
#>    gene            p_value weight_value
#>    <chr>             <dbl>        <dbl>
#>  1 leuC              0.523            1
#>  2 leuB              0.481            1
#>  3 leuA              0.508            1
#>  4 leuL              1                1
#>  5 leuO              0.802            1
#>  6 yabR              1                1
#>  7 ilvI              0.732            1
#>  8 ilvH              0.729            1
#>  9 BW25113_RS00395   1                1
#> 10 cra               0.571            1
#> # ℹ 21 more rows
```

Using a simple “Bonferroni correction” with $\alpha=0.05$ for multiple
testing problem ConNIS declared the follwing 13 genes as essential:

``` r
results_ConNIS %>% filter(p_value <= 0.05/nrow(truncated_ecoli))
#> # A tibble: 13 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 ftsI  2.65e-11            1
#>  2 murE  1.66e-14            1
#>  3 murF  5.23e- 7            1
#>  4 mraY  4.75e- 9            1
#>  5 murD  8.54e- 8            1
#>  6 ftsW  2.81e- 5            1
#>  7 murG  4.95e- 9            1
#>  8 murC  1.73e-14            1
#>  9 ftsQ  2.91e- 6            1
#> 10 ftsA  1.58e- 5            1
#> 11 ftsZ  3.95e- 9            1
#> 12 lpxC  7.81e- 9            1
#> 13 secA  6.30e- 6            1
```

Next, we re-run ConNIS with a smaller weight and apply again a the
“Bonferroni” method.

``` r
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 0.2)

results_ConNIS %>% filter(p_value <= 0.05/nrow(truncated_ecoli))
#> # A tibble: 8 × 3
#>   gene   p_value weight_value
#>   <chr>    <dbl>        <dbl>
#> 1 ftsI  0.000339          0.2
#> 2 murE  0.000672          0.2
#> 3 mraY  0.000924          0.2
#> 4 murG  0.000937          0.2
#> 5 murC  0.000678          0.2
#> 6 ftsQ  0.00120           0.2
#> 7 ftsZ  0.000869          0.2
#> 8 lpxC  0.00109           0.2
```

Only 8 genes are declared essential since smaller weights will make it
harder to label a gene (by chance) as ‘essential’.

Next, we give an example how to select a weight by the instability
approach. We will use the parallel version of the function using
`mclapply`. We use different weights. NOTE: This might take a bit since
we use 200 subsamples. You can reduce the number of subsamples.

``` r

# set weight
weights <- seq(0.1, 1, 0.1)

out <- 
  instabilities(
  method="ConNIS", 
  sig.level = 0.05, 
  p.adjust.mehtod = "bonferroni", 
  ins.positions = truncated_is_pos, 
  gene.names = truncated_ecoli$gene, 
  gene.starts = truncated_ecoli$start, 
  gene.stops = truncated_ecoli$end, 
  genome.length = max(truncated_ecoli$end), 
  weights = weights, 
  m = 30, 
  d = 0.5, 
  use.parallelization = T, 
  parallelization.type = "mclapply", 
  numCores = detectCores()-1, 
  seed = 1, 
  set.rng = "L'Ecuyer-CMRG")

out
#> # A tibble: 10 × 2
#>    weight_value instability
#>           <dbl>       <dbl>
#>  1          0.1       0.157
#>  2          0.2       0.157
#>  3          0.3       0.157
#>  4          0.4       0.147
#>  5          0.5       0.147
#>  6          0.6       0.195
#>  7          0.7       0.195
#>  8          0.8       0.195
#>  9          0.9       0.195
#> 10          1         0.195
```

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
