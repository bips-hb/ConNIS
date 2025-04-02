
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ConNIS <img src="./man/figures/logo.svg" alt="ConNIS" align="right" width="120"/>

<!-- badges: start -->
<!-- badges: end -->

The `ConNIS` package provides the implementation of the *Consecutive
Non-Insertion Sites* (ConNIS) method (Hanke et al., 2025) which
determines putative essential genes based on the probabilities to
observe sequences of non-insertions within genes by chance. The method
allows to set a weight value $w$ to adjust for non-uniformly distributed
gene-wise insertion sites across the genome. For a data driven selection
of $w$ an labeling instability approach is provided.

In addition, the following methods are implemented in the package:

- the *Binomial* approach of the `TSAS 2.0` package (Burger et
  al., 2017) with an additional weigh parameter

- the *Tn5Gaps* method of the `TRANSIT` package (DeJesus et al, 2015)
  with an additional weigh parameter

- the determination of essential genes by applying the *Geometric*
  distribution with an additional weight parameter

Each of this method can be adjusted by $w$, too.

## Installation

You can install the development version of ConNIS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("moritz-hanke/ConNIS")
```

## Example

We use a truncated real world dataset and a list of truncated insertion
sites as an toy example for applying ConNIS.

``` r

library(ConNIS)

# Use the E. coli BW 25113 dataset but only the first 30 genes
truncated_ecoli <- ecoli_bw25113[1:30,]

# load the insertion sites by Goodall et al., 2018, but omit all insertion sites 
# that do not lie within the genomic region of truncated_ecoli
truncated_is_pos <- sort(is_pos[is_pos <= max(truncated_ecoli$end) &
                           is_pos >= min(truncated_ecoli$start)])

# run ConNIS with weight = 1
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 1)

results_ConNIS
#> # A tibble: 30 × 3
#>    gene  p_value weight_value
#>    <chr>   <dbl>        <dbl>
#>  1 thrL  0.446              1
#>  2 thrA  0.0971             1
#>  3 thrB  0.166              1
#>  4 thrC  0.00429            1
#>  5 yaaX  0.326              1
#>  6 yaaA  0.0142             1
#>  7 yaaJ  0.0130             1
#>  8 talB  0.0896             1
#>  9 mog   0.0138             1
#> 10 satP  0.0265             1
#> # ℹ 20 more rows
```

Using the “Bonferroni-Holm” correction with $\alpha=0.05$ for multiple
testing adjustment ConNIS declared the follwing 12 genes as essential:

``` r
results_ConNIS %>% filter(p.adjust(p_value, "BH") <= 0.05)
#> # A tibble: 12 × 3
#>    gene   p_value weight_value
#>    <chr>    <dbl>        <dbl>
#>  1 thrC  4.29e- 3            1
#>  2 yaaA  1.42e- 2            1
#>  3 yaaJ  1.30e- 2            1
#>  4 mog   1.38e- 2            1
#>  5 dnaK  6.73e-14            1
#>  6 rpsT  5.13e- 9            1
#>  7 ribF  6.29e-35            1
#>  8 ileS  4.46e-24            1
#>  9 lspA  1.20e-17            1
#> 10 fkpB  2.78e- 3            1
#> 11 ispH  1.58e- 7            1
#> 12 dapB  1.54e-30            1
```

Next, we re-run ConNIS with a smaller weight value and apply again a the
“Bonferroni” method.

``` r
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 0.2)

results_ConNIS %>% filter(p.adjust(p_value, "BH") <= 0.05)
#> # A tibble: 6 × 3
#>   gene   p_value weight_value
#>   <chr>    <dbl>        <dbl>
#> 1 dnaK  2.08e- 3          0.2
#> 2 rpsT  3.80e- 3          0.2
#> 3 ribF  3.08e-11          0.2
#> 4 ileS  2.45e- 5          0.2
#> 5 lspA  8.21e- 6          0.2
#> 6 dapB  1.09e- 7          0.2
```

6 genes are declared essential since smaller weights will make it harder
to label a gene (by chance) as ‘essential’.

Next, we give an example for selecting `weight` our of five different
values by the instability approach. We will use the function’s build-in
parallelization option (`parallelization.type = "mclapply"`). NOTE: Only
10 subsamples are drawn for demonstration purpose. For real world
applications we suggest $m \approx  500$.

``` r

# set weight
weights <- seq(0.2, 1, 0.2)

instabilities_connis <- 
  instabilities(
  method="ConNIS", 
  sig.level = 0.05, 
  p.adjust.mehtod = "BH", 
  ins.positions = truncated_is_pos, 
  gene.names = truncated_ecoli$gene, 
  gene.starts = truncated_ecoli$start, 
  gene.stops = truncated_ecoli$end, 
  genome.length = max(truncated_ecoli$end), 
  weights = weights, 
  m = 10, 
  d = 0.5, 
  use.parallelization = T, 
  parallelization.type = "mclapply", 
  numCores = detectCores()-1, 
  seed = 1, 
  set.rng = "L'Ecuyer-CMRG")

instabilities_connis
#> # A tibble: 5 × 2
#>   weight_value instability
#>          <dbl>       <dbl>
#> 1          0.2       0.148
#> 2          0.4       0.16 
#> 3          0.6       0.138
#> 4          0.8       0.149
#> 5          1         0.145
```

Applying ConNIS again with the weight value that had the minimal
instability, i.e. `weight = 0.6`, 7 genes are declared essential:

``` r
results_ConNIS <- 
  ConNIS(ins.positions = truncated_is_pos, 
       gene.names = truncated_ecoli$gene, 
       gene.starts = truncated_ecoli$start, 
       gene.stops = truncated_ecoli$end, 
       genome.length = max(truncated_ecoli$end), 
       weight = 0.6)

results_ConNIS %>% filter(p.adjust(p_value, "BH") <= 0.05)
#> # A tibble: 7 × 3
#>   gene   p_value weight_value
#>   <chr>    <dbl>        <dbl>
#> 1 dnaK  8.12e- 9          0.6
#> 2 rpsT  3.34e- 7          0.6
#> 3 ribF  6.99e-24          0.6
#> 4 ileS  1.20e-14          0.6
#> 5 lspA  4.16e-12          0.6
#> 6 ispH  4.19e- 5          0.6
#> 7 dapB  2.24e-21          0.6
```

## References

- `Burger, B. T., Imam, S., Scarborough, M. J., Noguera, D. R. & Donohue, T. J. Combining genome-scale experimental and computational methods to identify essential genes in rhodobacter sphaeroides. mSystems 2 (2017). URL http://dx. doi.org/10.1128/msystems.00015-17`
- `DeJesus, M. A., Ambadipudi, C., Baker, R., Sassetti, C. & Ioerger, T. R. Transit - a software tool for himar1 tnseq analysis. PLOS Computational Biology 11, e1004401 (2015). URL http://dx.doi.org/10.1371/journal.pcbi.1004401`
- `Goodall, E. C. A. et al. The essential genome of escherichia coli k-12. mBio 9 (2018). URL http://dx.doi.org/10.1128/mBio.02096-17.`
