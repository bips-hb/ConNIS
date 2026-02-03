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
                   weight = 1) {
  
  # Input checks
  if(length(gene.names) != length(gene.starts) || length(gene.names) != length(gene.stops)) {
    stop("Different lengths of gene.names, gene.starts and gene.stops")
  }
  
  gene.starts <- as.integer(gene.starts)
  gene.stops  <- as.integer(gene.stops)
  
  ins_sites <- sort(unique(as.integer(ins.positions)))
  M <- length(ins_sites)
  G <- length(gene.names)
  
  observed_genome_insertion_density <- M / genome.length
  
  # find slice indices per gene in O(G log M) 
  # For each gene interval [start, stop], compute indices of ins_sites within it:
  # left = first index with ins_sites >= start
  # right = last  index with ins_sites <= stop
  if(M > 0L) {
    left  <- findInterval(gene.starts - 1L, ins_sites) + 1L
    right <- findInterval(gene.stops,       ins_sites)
  } else {
    # No insertion sites in the entire genome
    left  <- rep.int(1L, G)
    right <- rep.int(0L, G)
  }
  
  # compute num.ins.per.gene with slice indices 
  if(is.null(num.ins.per.gene)) {
    num.ins.per.gene <- pmax(0L, right - left + 1L)
  } else {
    if(length(num.ins.per.gene) != G) {
      stop("Different lengths of gene.names, gene.starts, gene.stops and num.ins.per.gene")
    }
    num.ins.per.gene <- as.integer(num.ins.per.gene)
  }
  
  
  p_values <- numeric(G)
  
  # Helper: compute f(k) = 1 / choose(n-1, k-1) using gmp::chooseZ (big integer)
  f_at_k <- function(n, k) {
    # If k == 1, then choose(n-1, 0)=1 => f(k)=1
    if(k <= 1L) return(1)
    as.numeric(1 / chooseZ(n - 1L, k - 1L))
  }
  
  # Per-gene loop
  for (i in seq_len(G)) {
    
    start <- gene.starts[i]
    stop  <- gene.stops[i]
    n     <- as.integer(stop - start + 1L)        # gene length
    ins_n <- num.ins.per.gene[i]                  # observed insertions in gene
    
    # Expected number of insertion sites under non-essentiality
    expected_num_IS <- as.integer(floor(n * observed_genome_insertion_density))
    
    # k = number of empty urns
    k <- as.integer(n - expected_num_IS * weight)
    
    # The urn model requires 1 <= k <= n.
    if(k < 1L){
      k <- 1L
    }
    if(k > n){
      k <- n
    }
    
    # Compute max_gap using slice
    # use indices [left[i], right[i]] precomputed via findInterval().
    if(left[i] <= right[i]){
      ins_g <- ins_sites[left[i]:right[i]]
      gaps  <- diff(c(start - 1L, ins_g, stop + 1L))
    }else{
      # No insertions in this gene interval
      gaps <- n + 1L  
    }
    
    max_gap <- max(gaps)
    if(max_gap > n){
      max_gap <- n
    }  
    
    # --- Handle trivial case: insertion at every base of the gene (as in your original code) ---
    if(n == ins_n) {
      p_values[i] <- 1
      next
    }
    
    # If max_gap >= k, then p-value is just P(S = k) (tail is single point).
    # compute f(k) once via chooseZ (big integer safe), no full prob vector.
    if(max_gap >= k) {
      p_values[i] <- f_at_k(n, k)
      next
    }
    
    # Compute p-value
    # p_value = P(S >= max_gap) = sum_{s=max_gap}^{k} f(s, n, k)
    # If max_gap is small, compute complement: 1 - sum_{s=1}^{max_gap-1} f(s)
    # If max_gap is large, compute the tail directly via backward recurrence from f(k)
    
    if(max_gap <= 1000L) {
      
      # Complement approach: p = 1 - sum_{s=1}^{max_gap-1} f(s)
      # Start value:
      # f(1) = C(n-2, k-1) / C(n-1, k-1) = (n - k) / (n - 1)
      # (simple ratio, no gmp required)
      if(n <= 1L) {
        # degenerate gene length
        p_values[i] <- 1
        next
      }
      
      f_s <- (n - k) / (n - 1L)
      sum_small <- f_s
      
      # Recurrence:
      # f(s+1) = f(s) * (k - s) / (n - s - 1)
      # accumulate s=1..(max_gap-2) to cover up to (max_gap-1)
      if(max_gap > 2L) {
        for (s in 1L:(max_gap - 2L)) {
          f_s <- f_s * (k - s) / (n - s - 1L)
          sum_small <- sum_small + f_s
        }
      }
      
      p_values[i] <- 1 - sum_small
      
    }else{
      
      # Tail approach (reverse the problem) without generating all probabilities:
      # Compute f(k) = 1 / C(n-1, k-1) via chooseZ, then walk backwards:
      # f(s) = f(s+1) * (n - s - 1) / (k - s)
      # This sums only (k - max_gap + 1) terms.
      f_next <- f_at_k(n, k)  # f(k)
      tail_sum <- f_next
      
      # s = k-1 down to max_gap
      if(k - 1L >= max_gap) {
        for (s in (k - 1L):max_gap) {
          f_next <- f_next * (n - s - 1L) / (k - s)
          tail_sum <- tail_sum + f_next
        }
      }
      
      p_values[i] <- tail_sum
    }
  }
  
  tibble::tibble(
    gene = gene.names,
    p_value = p_values,
    weight_value = weight
  )
}

