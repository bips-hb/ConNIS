
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
truncated_ecoli <- ecoli_bw25113[501:1000,]

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
#>    gene            p_value weight_value
#>    <chr>             <dbl>        <dbl>
#>  1 tesA             0.422             1
#>  2 ybbA             0.372             1
#>  3 ybbP             0.134             1
#>  4 rhsD             0.183             1
#>  5 ybbC             0.403             1
#>  6 BW25113_RS02575  0.154             1
#>  7 BW25113_RS24960  0.0474            1
#>  8 BW25113_RS25970  0.584             1
#>  9 ylbG             0.465             1
#> 10 selU             0.0484            1
#> # ℹ 490 more rows
```

Using a simple “Bonferroni correction” with $\alpha=0.05$ for multiple
testing problem ConNIS declared the follwing 13 genes as essential:

``` r
results_ConNIS %>% filter(p_value <= 0.05/nrow(truncated_ecoli))
#> # A tibble: 29 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 lpxH  5.07e-17            1
#>  2 cysS  9.74e-34            1
#>  3 ybcJ  4.47e- 5            1
#>  4 folD  5.75e- 5            1
#>  5 lipA  8.47e-14            1
#>  6 mrdB  1.30e-26            1
#>  7 mrdA  6.66e-46            1
#>  8 nadD  1.06e-14            1
#>  9 holA  2.79e-24            1
#> 10 lptE  1.92e-14            1
#> # ℹ 19 more rows
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
#> # A tibble: 15 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 lpxH  3.84e- 6          0.2
#>  2 cysS  2.26e- 9          0.2
#>  3 mrdB  4.38e- 9          0.2
#>  4 mrdA  1.84e-12          0.2
#>  5 nadD  4.88e- 6          0.2
#>  6 holA  1.88e- 6          0.2
#>  7 lptE  5.94e- 6          0.2
#>  8 leuS  2.44e-18          0.2
#>  9 lnt   7.56e- 5          0.2
#> 10 lolA  5.37e- 6          0.2
#> 11 serS  3.71e- 6          0.2
#> 12 msbA  2.58e-12          0.2
#> 13 kdsB  3.60e- 6          0.2
#> 14 mukF  2.60e- 9          0.2
#> 15 asnS  2.19e- 9          0.2
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
#> 1          0.2      0.166 
#> 2          0.4      0.133 
#> 3          0.6      0.0952
#> 4          0.8      0.0867
#> 5          1        0.0898
```

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
