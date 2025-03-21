#' The Geometric analysis method
#'
#' Calculates the probability for each gene to contain a sequence of non-
#' insertion sites based on the geometric distribution
#' \deqn{ (1- \theta \cdot w) ^{l_j} (\theta \cdot w) .}
#' \eqn{l_j} is the length of the longest insertion-free sequence of gene
#' \eqn{j}. \eqn{\theta} is the genome-wide insertion density which can be
#' weighted by \eqn{0 < w \leq 1} to account for regions with low insertion
#' densities.
#'
#' @param ins.positions Numeric vector of position of observed insertions sites.
#' @param gene.names The names of the genes.
#' @param gene.starts Starting position within the genome of each gene.
#' @param gene.stops Ending position within the genome of each gene.
#' @param num.ins.per.gene Number of unique insertion sites within each gene.
#' @param genome.length Length of the genome.
#' @param weight A weight value for the genome-wide insertion density.
#'
#' @returns The p-values for each gene to observe its biggest gap.
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
#' Geometric(ins.positions = random_is, gene.names = genes, gene.starts = starts,
#' gene.stops = stops, genome.length = genome_length)
#'
#' @references
#' \insertRef{Burger2017}{ConNIS}
#'
#' @export

Geometric <- function(ins.positions,
                      gene.names,
                      gene.starts,
                      gene.stops,
                      num.ins.per.gene=NULL,
                      genome.length,
                      weight=1){

  if (is.null(num.ins.per.gene)) {
    if (!length(unique(
      c(
        length(gene.names),
        length(gene.starts),
        length(gene.stops)
      )
    )) == 1) {
      stop("Different lengths of gene.names, gene.starts, gene.stops and num.ins.per.gene")
    }

    num.ins.per.gene <- sapply(seq(gene.starts), function(start_i) {
      sum(ins.positions >= gene.starts[start_i] &
            ins.positions <= gene.stops[start_i])
    })
  } else {
    if (!length(unique(
      c(
        length(gene.names),
        length(gene.starts),
        length(gene.stops),
        length(num.ins.per.gene)
      )
    )) == 1) {
      stop("Different lengths of gene.names, gene.starts, gene.stops and num.ins.per.gene")
    }
  }

  ins_sites <- sort(ins.positions)

  observed_genome_insertion_densitiy <-
    length(ins_sites)/genome.length

  results_per_gene <- lapply(seq_along(gene.names), function(i){

    gene_i_start <-
      gene.starts[i]

    gene_i_stop <-
      gene.stops[i]

    gene_i_length <- gene_i_stop - gene_i_start + 1

    expected_num_IS_gene_i <-
      floor(gene_i_length * observed_genome_insertion_densitiy)

    gene_i_num_ins <- num.ins.per.gene[i]

    max_gap <- max(
      diff(
        unique(
          c(
            gene_i_start-1,
            ins_sites[ins_sites >= gene_i_start & ins_sites <= gene_i_stop],
            gene_i_stop+1
          )
        )
      )
    )
    max_gap <- min(gene_i_length, max_gap)

    p_value <- 1 - pgeom(max_gap, observed_genome_insertion_densitiy * weight)

    tibble(gene = gene.names[i],
           p_value = p_value,
           weight_value = weight)
  })

  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
}
