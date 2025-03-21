#' Instabilities for each weight based on subsamples
#'
#' 'instabilites' calculates the instabilities for list of weights for a given
#' methods. Basically, it is a wrapper around `subsampleResults()`.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @importFrom dplyr reframe
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom gmp chooseZ
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterSetRNGStream
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom Rdpack reprompt
#' @importFrom stats dbinom
#' @importFrom stats p.adjust
#' @importFrom stats pgeom
#' @importFrom utils menu
#' @importFrom tibble tibble
#'
#' @param method Method that will be applied to subsamples. Available are
#' `"Binomial"`, `"ConNIS"`, `"Geometric"` and `"Tn5Gaps"`.
#' @param sig.level Significance level that is applied for determining if a gene
#' is essentially
#' @param p.adjust.mehtod Adjustment method for the p.value for handeling
#' multiple testing. It uses the `stats` R function `p.adjust()`. Available values
#' are `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`,
#' `"fdr"` and `"none"`. See `?p.adjust` for further details.
#' @param ins.positions Numeric vector of position of observed insertions sites.
#' @param gene.names The names of the genes.
#' @param gene.starts Starting position within the genome of each gene.
#' @param gene.stops Ending position within the genome of each gene.
#' @param num.ins.per.gene Number of unique insertion sites within each gene.
#' @param genome.length Length of the genome.
#' @param weights A sequence of weight values that are applied to the
#' selected method.
#' @param m Number of sub samples.
#' @param d Proportion of the original IS used as sub samples.
#' @param use.parallelization Should the calculations be done in parallel?
#' @param parallelization.type Which method should be used for the
#' parallelization? Available are `"mclapply"` (fork mechanism) and
#' `"parLapply"` (creation of sockets) of the `"parallel"` package.
#' @param numCores Number of cores to juse for parallel calculation. Should not
#' exit the number of available (logical) cores of the system.
#' @param cluster.type If `"parLapply"` is used as parallelization.type, a
#' cluster type like `"PSOCK"` needs tp be specified. Check which types are
#' supported on your system. See als `?makeCluster` of the `parallel` package.
#' @param seed Seed for subsampling. If `NULL` no seed is set and results might
#' not be reproducable.
#' @param set.rng Which random number generator (RNG) should be used. If NULL, the
#' default RNG of the system is used. NOTE: if use.parallelizationization is `TRUE` the
#' `"L'Ecuyer-CMRG"` value for RNG should be used. For details see the manual of
#' the `parallel` package.
#' @param keep.new.RNG If the RNG is changed, should it kept for the rest of the
#' session? Default is `FALSE` (recommended).
#'
#' @returns A `tibble` with the instability for each weight.
#'
#' @examples
#' # generate random insertion sites
#' set.seed(1)
#' random_is <- sort(sample(1:10000, 2000))
#'
#' # generate dummy gene names
#' genes <- paste("gene_", 1:30, sep="")
#'
#' # generate radom start points and stop points of the dummy genes
#' set.seed(2)
#' x <- sort(sample(1:10000, 60))
#' starts <- x[seq(1,60, 2)]
#' stops <- x[seq(2,60, 2)]
#'
#' genome_length <- 10000
#'
#' # run instability approach with ConNIS (sequential)
#' instabilities(method="ConNIS",
#'                   sig.level = 0.05,
#'                   p.adjust.mehtod = "bonferroni",
#'                   ins.positions = random_is,
#'                   gene.names = genes,
#'                   gene.starts = starts,
#'                   gene.stops = stops,
#'                   genome.length = genome_length,
#'                   weights = c(0.3, 0.6, 0.9),
#'                   m = 10,
#'                   d = 0.5,
#'                   use.parallelization = FALSE,
#'                   set.rng = "L'Ecuyer-CMRG",
#'                   seed = 1)
#'
#' # Example with parallelization (mclapply()); use detectCores()-1 to set
#' # workers
#' \dontrun{
#' instabilities(method="ConNIS",
#'                   sig.level = 0.05,
#'                   p.adjust.mehtod = "bonferroni",
#'                   ins.positions = random_is,
#'                   gene.names = genes,
#'                   gene.starts = starts,
#'                   gene.stops = stops,
#'                   genome.length = genome_length,
#'                   weights = c(0.3, 0.6, 0.9),
#'                   m = 100,
#'                   d = 0.5,
#'                   use.parallelization = T,
#'                   parallelization.type = "mclapply",
#'                   set.rng = "L'Ecuyer-CMRG",
#'                   numCores = max(1,detectCores()-1),
#'                   seed = 1)
#' }
#' @export

instabilities <- function(
    method="ConNIS",
    sig.level = 0.05,
    p.adjust.mehtod = "bonferroni",
    ins.positions,
    gene.names,
    gene.starts,
    gene.stops,
    num.ins.per.gene = NULL,
    genome.length,
    weights = 1,
    m = 100,
    d = 0.5,
    use.parallelization = FALSE,
    parallelization.type = "mclapply",
    numCores = NULL,
    cluster.type = NULL,
    seed = NULL,
    set.rng = NULL,
    keep.new.RNG=FALSE){

  list_of_instabilities <- lapply(weights, function(w){

    sub_results <- subsampleResults(
      method = method,
      ins.positions = ins.positions,
      gene.names = gene.names,
      gene.starts = gene.starts,
      gene.stops = gene.stops,
      num.ins.per.gene = num.ins.per.gene,
      genome.length = genome.length,
      weights = w,
      m = m,
      d = d,
      use.parallelization = use.parallelization,
      parallelization.type = parallelization.type,
      numCores = numCores,
      cluster.type = cluster.type,
      seed = seed,
      set.rng = set.rng,
      keep.new.RNG = keep.new.RNG
    )

    do.call(rbind, sub_results) %>%
      group_by(.data$subsample_number) %>%
      mutate(p_value_adjusted = p.adjust(.data$p_value, method = p.adjust.mehtod)) %>%
      mutate(ess = .data$p_value_adjusted <= sig.level) %>%
      group_by(.data$gene) %>%
      reframe(mean_ess = sum(.data$ess)/length(unique(.data$subsample_number))) %>%
      reframe(var_ess = .data$mean_ess*(1-.data$mean_ess)) %>%
      reframe(weight_value =w,
              instability = sum(.data$var_ess)/sum(.data$var_ess>0))

  })

  do.call(rbind, list_of_instabilities)

}


