#' Number of possible sequences of empty urns
#'
#' `prob_seq_misses` calculates the probability to observe *s* adjacent empty
#' urns if *n* urns are placed in a row of which *k* urns are randomly
#' filled:
#' \deqn{
#' \frac{\binom{n-s-2}{k-s}}{\binom{n-1}{k-1}}
#' }
#'
#' @import gmp
#'
#' @param n Length of sequence.
#' @param k Number total misses.
#'
#' @returns A list of length k with number of possibilities to observe a sequence of misses of length 1, 2, ..., k.
#'
#' @examples
#' prob_seq_misses(23,8)
#'
#' @export

prob_seq_misses <- function(n, k) {
  counts <- freq_seq_misses(n, k)

  freqs <- Reduce("+", counts)
  all_probs <-
    sapply(seq_along(counts), function(i) {
      as.numeric(counts[[i]] / freqs)
    })
  all_probs
}
