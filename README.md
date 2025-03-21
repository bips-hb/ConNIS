
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
truncated_ecoli <- ecoli_bw25113[50:100,]

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
#> # A tibble: 51 × 3
#>    gene             p_value weight_value
#>    <chr>              <dbl>        <dbl>
#>  1 pdxA            1.77e- 1            1
#>  2 surA            1.91e- 1            1
#>  3 lptD            4.88e-35            1
#>  4 djlA            2.63e- 1            1
#>  5 yabP            4.78e- 1            1
#>  6 yabQ            5.51e- 1            1
#>  7 BW25113_RS25400 7.75e- 1            1
#>  8 rluA            3.55e- 1            1
#>  9 rapA            1.76e- 1            1
#> 10 polB            1.22e- 1            1
#> # ℹ 41 more rows
```

Using a simple “Bonferroni correction” with $\alpha=0.05$ for multiple
testing problem ConNIS declared the follwing 13 genes as essential:

``` r
results_ConNIS %>% filter(p_value <= 0.05/nrow(truncated_ecoli))
#> # A tibble: 14 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 lptD  4.88e-35            1
#>  2 ftsI  2.84e-28            1
#>  3 murE  3.78e-32            1
#>  4 murF  1.73e-18            1
#>  5 mraY  1.85e-22            1
#>  6 murD  5.07e-21            1
#>  7 ftsW  3.06e-12            1
#>  8 murG  2.09e-22            1
#>  9 murC  4.72e-30            1
#> 10 ftsQ  1.91e-17            1
#> 11 ftsA  3.91e-14            1
#> 12 ftsZ  9.25e-25            1
#> 13 lpxC  8.32e-20            1
#> 14 secA  4.34e-14            1
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
#> # A tibble: 12 × 3
#>    gene        p_value weight_value
#>    <chr>         <dbl>        <dbl>
#>  1 lptD  0.000000267            0.2
#>  2 ftsI  0.00000595             0.2
#>  3 murE  0.00000000183          0.2
#>  4 murF  0.0000203              0.2
#>  5 mraY  0.00000171             0.2
#>  6 murD  0.00000534             0.2
#>  7 murG  0.00000176             0.2
#>  8 murC  0.00000000187          0.2
#>  9 ftsQ  0.00000291             0.2
#> 10 ftsA  0.000255               0.2
#> 11 ftsZ  0.00000151             0.2
#> 12 lpxC  0.00000238             0.2
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
#> 1          0.2       0.195
#> 2          0.4       0.201
#> 3          0.6       0.154
#> 4          0.8       0.154
#> 5          1         0.154
```

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
