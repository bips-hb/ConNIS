#' Function to calculate the number of possibles to observe sequences of misses
#'
#' @import gmp
#' @param n Length of sequence.
#' @param k Number total misses.
#' @returns A list of length k with number of possibilities to observe a sequence of misses of length 1, 2, ..., k.
#' @examples
#' freq_seq_misses(23,8)
#' @export

freq_seq_misses <- function(n, k) {
  lapply(1:k, function(s) {
    (n - k + 1) * gmp::chooseZ(n - s - 1, k - s)
  })
}
