
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

We use a truncated real world dataset and a list of truncated insertion
sites as an example for applying ConNIS.

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
truncated_ecoli <- ecoli_bw25113[1:500,]

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
#> # A tibble: 500 × 3
#>    gene  p_value weight_value
#>    <chr>   <dbl>        <dbl>
#>  1 thrL  0.446              1
#>  2 thrA  0.0568             1
#>  3 thrB  0.106              1
#>  4 thrC  0.00129            1
#>  5 yaaX  0.260              1
#>  6 yaaA  0.00376            1
#>  7 yaaJ  0.00454            1
#>  8 talB  0.0505             1
#>  9 mog   0.00578            1
#> 10 satP  0.0127             1
#> # ℹ 490 more rows
```

Using a simple “Bonferroni correction” with $\alpha=0.05$ for multiple
testing problem ConNIS declared the follwing 13 genes as essential:

``` r
results_ConNIS %>% filter(p_value <= 0.05/nrow(truncated_ecoli))
#> # A tibble: 73 × 3
#>    gene    p_value weight_value
#>    <chr>     <dbl>        <dbl>
#>  1 dnaK  4.71e- 17            1
#>  2 rpsT  9.91e- 11            1
#>  3 ribF  1.01e- 41            1
#>  4 ileS  8.54e- 30            1
#>  5 lspA  4.60e- 21            1
#>  6 ispH  9.50e-  9            1
#>  7 dapB  2.69e- 37            1
#>  8 folA  6.27e- 21            1
#>  9 lptD  9.96e-101            1
#> 10 ftsL  3.61e-  6            1
#> # ℹ 63 more rows
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
#> # A tibble: 45 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 ribF  3.08e-11          0.2
#>  2 ileS  1.00e- 6          0.2
#>  3 lspA  8.21e- 6          0.2
#>  4 dapB  2.66e-10          0.2
#>  5 folA  8.74e- 6          0.2
#>  6 lptD  3.98e-19          0.2
#>  7 ftsI  6.01e-15          0.2
#>  8 murE  3.18e-19          0.2
#>  9 murF  3.20e-10          0.2
#> 10 mraY  8.17e-14          0.2
#> # ℹ 35 more rows
```

Only 8 genes are declared essential since smaller weights will make it
harder to label a gene (by chance) as ‘essential’.

Next, we give an example how to select a weight by the instability
approach. We will use the parallel version of the function using
`mclapply`. We use different weights. NOTE: We use only 30 subsamples
for demonstration purpose. For real world applications we suggest
$m \approx  500$ and parallel processing uization of the .

``` r

# set weight
weights <- seq(0.2, 1, 0.2)

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
#> # A tibble: 5 × 2
#>   weight_value instability
#>          <dbl>       <dbl>
#> 1          0.2      0.178 
#> 2          0.4      0.129 
#> 3          0.6      0.116 
#> 4          0.8      0.0941
#> 5          1        0.0880
```

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
