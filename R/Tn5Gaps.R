#' The Tn5Gaps analysis method of TRANSIT
#'
#' Implementation of the TntGaps method by
#' \insertCite{DeJesus2015;textual}{ConNIS} which uses the Gumbel distribution.
#' See https://transit.readthedocs.io/en/v3.2.8/method_tn5gaps.html for further
#' details and https://github.com/ioerger/transit for the original Python
#' implementation.
#' Here, the method is extended by allowing to set a weight value
#' \eqn{0 < w \leq 1} to account for genomic regions with an insertion density
#' below the genome-wide average.
#'
#'
#' @param ins.positions Numeric vector of position of observed insertions sites.
#' @param gene.names The names of the genes.
#' @param gene.starts Starting position within the genome of each gene.
#' @param gene.stops Ending position within the genome of each gene.
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
#' Tn5Gaps(ins.positions = random_is, gene.names = genes, gene.starts = starts,
#' gene.stops = stops, genome.length = genome_length)
#'
#' @references
#' \insertRef{DeJesus2015}{ConNIS}
#'
#' @export

Tn5Gaps <- function(ins.positions,
                    gene.names,
                    gene.starts,
                    gene.stops,
                    genome.length,
                    weighting=1){

  if(!length(unique(
    c(length(gene.names),
      length(gene.starts),
      length(gene.stops))))==1){
    stop("Different lengths of gene.names, gene.starts and gene.stops")
  }

  ins_sites <-
    sort(unique(
      c(1,
        unique(ins.positions),
        genome.length)))

  pins <- length(unique(ins.positions))/genome.length * weighting

  pnon <- 1 - pins

  results_per_gene <- lapply(seq_along(gene.names), function(i){

    gene_i_start <-
      gene.starts[i]

    gene_i_stop <-
      gene.stops[i]

    gene_i_length <- gene_i_stop - gene_i_start + 1

    lower_ins_index <-
      max(
        which(ins_sites <= gene_i_start)
      )

    upper_ins_index <-
      min(
        which(ins_sites >= gene_i_stop)
      )
    if(upper_ins_index == Inf){
      upper_ins_index <- lower_ins_index
    }

    gaps_gene_i <- diff(ins_sites[lower_ins_index:upper_ins_index])

    max_gap_gene_i <- unique(max(gaps_gene_i))

    ins_sites_overlapping_gene_i <-
      ins_sites[lower_ins_index:upper_ins_index]

    overlap_sizes_gene_i <- sapply(1:(length(ins_sites_overlapping_gene_i)-1), function(j){
      min(ins_sites_overlapping_gene_i[j+1],
          gene_i_stop
      ) -
        max(ins_sites_overlapping_gene_i[j], gene_i_start) -1
    })

    max_overlap_gene_i <- unique(max(overlap_sizes_gene_i))

    gap_size_for_p_value <-
      gaps_gene_i[overlap_sizes_gene_i %in% max_overlap_gene_i][1]

    # pvalue from tn5gaps.pv of transit software
    # https://github.com/mad-lab/transit/blob/master/src/pytransit/analysis/tn5gaps.py
    p_value <-
      1 - pgumbel(
        k = gap_size_for_p_value,
        mu = log(length(ins_sites)* pins, 1/pnon),
        s = 1/log(1/pnon))


    tibble(gene = gene.names[i],
           p_value = p_value,
           weight_value = weighting)

  })

  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
}
