#' The ConNIS analysis method
#'
#' Calculates the probability for each gene of its biggest observed insertion
#' free gap. For a gene \eqn{j} ConNIS is defined as
#' \deqn{1-\sum_{i=1}^{l_j-1} \frac{\binom{b_j-i-1}{b_j - \hat{h}_j w-i}}{\binom{b_j-1}{b_j - \hat{h}_j w -1}}}
#' with \eqn{b_j} being the length of gene \eqn{j} and \eqn{l_j} being its
#' longest sequences without insertion sites. \eqn{ \hat{h}_j = b_j \cdot \theta}
#' are the estimated number of insertions under the assumption of
#' non-essentiality based on the genome-wide insertion densitiy \eqn{\theta}.
#' \eqn{w} we is a weight value to adjust for regions with an insertion density
#' below \eqn{\theta}.
#'
#' @importFrom Rdpack reprompt
#' @importFrom tibble tibble
#'
#' @param ins.positions Numeric vector of position of observed insertions sites.
#' @param gene.names The names of the genes.
#' @param gene.starts Starting position within the genome of each gene.
#' @param gene.stops Ending position within the genome of each gene.
#' @param num.ins.per.gene Number of unique insertion sites within each gene.
#' @param genome.length Length of the genome.
#' @param weighting A weighting value for the genome-wide insertion density.
#'
#' @returns The p-values for each gene to observe its biggest gap.
#'
#' @examples
#' # generate random insertion sites
#' set.seed(1)
#' random_is <- sort(sample(1:10000, 2000))
#'
#' # generate dummy gene names
#' genes <- paste("gene_", 1:30)
#'
#' # generate radom start points and stop points of the dummy genes
#' set.seed(2)
#' x <- sort(sample(1:10000, 60))
#' starts <- x[seq(1,60, 2)]
#' stops <- x[seq(2,60, 2)]
#'
#' genome_length <- 10000
#'
#' ConNIS(ins.positions = random_is, gene.names = genes, gene.starts = starts,
#' gene.stops = stops, genome.length = genome_length)
#'
#' @export

ConNIS <- function(ins.positions,
                   gene.names,
                   gene.starts,
                   gene.stops,
                   num.ins.per.gene = NULL,
                   genome.length,
                   weighting = 1) {
  if (is.null(num.ins.per.gene)) {
    if (!length(unique(
      c(
        length(gene.names),
        length(gene.starts),
        length(gene.stops)
      )
    )) == 1) {
      stop("Different lengths of gene.names, gene.starts and gene.stops")
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
    length(ins_sites) / genome.length

  results_per_gene <- lapply(seq_along(gene.names), function(i) {
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
            gene_i_start - 1,
            ins_sites[ins_sites >= gene_i_start & ins_sites <= gene_i_stop],
            gene_i_stop + 1
          )
        )
      )
    )
    max_gap <- min(gene_i_length, max_gap)

    # applying ConNIS for gaps <= 1000 can be speed up by reversing the problem

    if (max_gap > 1000) {
      probs <-
        prob_seq_misses(
          gene_i_length,
          gene_i_length -
            expected_num_IS_gene_i * weighting
        )
      p_value <-
        sum(probs[max_gap:(gene_i_length -
          expected_num_IS_gene_i * weighting)])

      if (gene_i_length == gene_i_num_ins) {
        p_value <- 1
      }
      if (max_gap >= gene_i_length - expected_num_IS_gene_i * weighting) {
        p_value <- probs[length(probs)]
      }
    } else {
      if (gene_i_length == gene_i_num_ins) {
        p_value <- 1
      }
      if (max_gap >= gene_i_length - expected_num_IS_gene_i * weighting) {
        p_value <-
          as.numeric(chooseZ(
            gene_i_length - (gene_i_length - expected_num_IS_gene_i * weighting) - 1,
            gene_i_length - expected_num_IS_gene_i * weighting -
              (gene_i_length - expected_num_IS_gene_i * weighting)
          ) /
            chooseZ(
              gene_i_length - 1,
              gene_i_length - expected_num_IS_gene_i * weighting - 1
            ))
      } else {
        probs <-
          sapply(1:(max_gap - 1), function(s) {
            as.numeric(chooseZ(
              gene_i_length - s - 1,
              gene_i_length - expected_num_IS_gene_i * weighting - s
            ) /
              chooseZ(
                gene_i_length - 1,
                gene_i_length - expected_num_IS_gene_i * weighting - 1
              ))
          })

        p_value <- 1 - sum(probs)
      }
    }

    tibble(
      gene = gene.names[i],
      p_value = p_value,
      weight_value = weighting
    )
  })

  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
}

