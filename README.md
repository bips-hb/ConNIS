
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
truncated_ecoli <- ecoli_bw25113[51:150,]

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
#> # A tibble: 100 × 3
#>    gene             p_value weight_value
#>    <chr>              <dbl>        <dbl>
#>  1 surA            7.69e- 2            1
#>  2 lptD            7.68e-58            1
#>  3 djlA            1.22e- 1            1
#>  4 yabP            3.06e- 1            1
#>  5 yabQ            3.02e- 1            1
#>  6 BW25113_RS25400 6.00e- 1            1
#>  7 rluA            1.54e- 1            1
#>  8 rapA            6.19e- 2            1
#>  9 polB            3.41e- 2            1
#> 10 araD            8.98e- 2            1
#> # ℹ 90 more rows
```

Using a simple “Bonferroni correction” with $\alpha=0.05$ for multiple
testing problem ConNIS declared the follwing 13 genes as essential:

``` r
results_ConNIS %>% filter(p_value <= 0.05/nrow(truncated_ecoli))
#> # A tibble: 21 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 lptD  7.68e-58            1
#>  2 ftsL  4.29e- 4            1
#>  3 ftsI  3.31e-47            1
#>  4 murE  9.90e-46            1
#>  5 murF  9.76e-31            1
#>  6 mraY  3.15e-32            1
#>  7 murD  6.91e-36            1
#>  8 ftsW  1.37e-20            1
#>  9 murG  3.83e-32            1
#> 10 murC  1.17e-45            1
#> # ℹ 11 more rows
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
#> # A tibble: 16 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 lptD  2.58e-12          0.2
#>  2 ftsI  1.67e- 9          0.2
#>  3 murE  4.93e-12          0.2
#>  4 murF  5.23e- 7          0.2
#>  5 mraY  4.75e- 9          0.2
#>  6 murD  8.54e- 8          0.2
#>  7 ftsW  2.81e- 5          0.2
#>  8 murG  4.95e- 9          0.2
#>  9 murC  5.09e-12          0.2
#> 10 ftsQ  1.05e- 8          0.2
#> 11 ftsA  1.58e- 5          0.2
#> 12 ftsZ  3.95e- 9          0.2
#> 13 lpxC  7.81e- 9          0.2
#> 14 secA  2.39e- 5          0.2
#> 15 lpd   5.86e-12          0.2
#> 16 can   4.57e- 6          0.2
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
  m = 100, 
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
#> 1          0.2      0.173 
#> 2          0.4      0.133 
#> 3          0.6      0.117 
#> 4          0.8      0.0890
#> 5          1        0.0800
```

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
