#' Probabilities of observing a sequence of empty urns of a given length
#'
#' `prob_seq_misses` calculates the probabilities to observe
#' \eqn{s = 1, \dots, k} adjacent empty urns if \eqn{n} urns are placed in a row
#' of which \eqn{k} urns are empty and randomly placed:
#' \deqn{
#' \frac{\binom{n-s-2}{k-s}}{\binom{n-1}{k-1}}
#' }
#'
#' @import gmp
#'
#' @param n Number of urns.
#' @param k Number empty urns.
#' @param s Optional. Calculate probabilities only for a give sequence of length
#' `s` of empty urns. Must not be larger than `k` or smaller than `1`. If `NULL`
#' all probabilities of lengths \eqn{s=1,2, \dots, k} are calculated.
#'
#' @returns A vector of length 'k' with number of possibilities to observe a
#' sequence of empty urns. First list element gives the probability for
#' \eqn{s=1}, second list element for \eqn{2} and so on. If `s` is not NULL only
#' the corresponding probability is returned.
#'
#' @examples
#' prob_seq_misses(23,8)
#'
#' @export

prob_seq_misses <- function(n, k, s=NULL) {

  if(!is.null(s)){

    if(s > k | s < 1){

      stop("s needs to be >0 and <=k.")

    }else{

      counts <- freq_seq_misses(n, k)

      freqs <- Reduce("+", counts)

      as.numeric(counts[[s]] / freqs)

    }

  }else{

    counts <- freq_seq_misses(n, k)

    freqs <- Reduce("+", counts)

    sapply(seq_along(counts), function(i) {
      as.numeric(counts[[i]] / freqs)
    })

  }

}
