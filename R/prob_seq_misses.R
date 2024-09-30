#' Function to calculate the probabilities to observe sequences of misses
#'
#' @import gmp
#' @param n Length of sequence.
#' @param k Number total misses.
#' @returns A vector of length k with probabilities to observe a sequence of misses of length 1, 2, ..., k.
#' @examples
#' prob_seq_misses(23,8)
#'
prob_seq_misses <- function(n, k) {
  counts <- freq_seq_misses(n, k)

  freqs <- Reduce("+", counts)
  all_probs <-
    sapply(seq_along(counts), function(i) {
      as.numeric(counts[[i]] / freqs)
    })
  all_probs
}
