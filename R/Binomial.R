#' Calculate the probability for each gene to contain a number of observed IS
#' based on the binomial distribution
#'
#' The method is an implementation of the Tn-seq analysis software (TSAS 2.0)
#' package \insertRef{Burger2017}{ConNIS} (see https://github.com/srimam/TSAS).
#'
#' @importFrom stats dbinom
#' @param ins.positions The observed unique insertion sites.
#' @param gene.names The names of the genes.
#' @param gene.starts Starting position within the genome of each gene.
#' @param gene.stops Ending position within the genome of each gene.
#' @param num.ins.per.gene Number of unique insertion sites within each gene.
#' @param genome.length Length of the genome.
#' @param weighting A weighting value for the genome-wide insertion density.
#' @returns The p-values for each gene to observe its biggest gap.
#' @examples
#' set.seed(1234)
#' random_is <- sort(sample(1:10000, 2000))
#' genes <- paste("gene_", 1:30)
#' set.seed(5678)
#' x <- sort(sample(1:10000, 60))
#' starts <- x[seq(1,60, 2)]
#' stops <- x[seq(2,60, 2)]
#' genome_length <- 10000
#' Binomial(ins.positions = random_is, gene.names = genes, gene.starts = starts,
#' gene.stops = stops, genome.length = genome_length)
#' @export

Binomial <- function(ins.positions,
                     gene.names,
                     gene.starts,
                     gene.stops,
                     num.ins.per.gene=NULL,
                     genome.length,
                     weighting=1){

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
            ins.positions <= gene.stops[start_i]) /
        (gene.stops[start_i] - gene.starts[start_i] + 1)
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
    gene_i_num_ins <- num.ins.per.gene[i]

    p_value <- sum(dbinom(0:gene_i_num_ins,
                          gene_i_length,
                          observed_genome_insertion_densitiy * weighting))

    if(is.na(p_value)){
      p_value <- 1
    }

    tibble(gene = gene.names[i],
           p_value = p_value,
           weight_insdens = weighting)


  })

  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene


}

