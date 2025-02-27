#' Number of possible sequences of empty urns
#'
#' `freq_seq_misses` calculates the number of possibles to observe *s* adjacent
#' empty urns if *n* urns are placed in a row of which *k* urns are randomly
#' filled:
#' \deqn{
#' (n-s-1)\binom{n-s-2}{k-s}
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
#' freq_seq_misses(23,8)
#'
#' @export

freq_seq_misses <- function(n, k) {
  lapply(1:k, function(s) {
    (n - k + 1) * gmp::chooseZ(n - s - 1, k - s)
  })
}
