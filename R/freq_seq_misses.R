#' Number of possibilities of observing a sequence of empty urns of a given
#' length
#'
#' `freq_seq_misses` calculates the number of all possibles to observe
#' \eqn{s = 1, \dots, k} adjacent empty urns if \eqn{n} urns are placed in a row
#' of which \eqn{k} urns are empty and randomly placed:
#' \deqn{
#' (n-s-1)\binom{n-s-2}{k-s}
#' }
#' The function uses the `chooseZ` function of the `gmp` package.
#'
#' @param n Number of urns.
#' @param k Number of empty urns.
#' @param s Optional. Calculate possibilities only for a give sequence of length
#' `s` of empty urns. Must not be larger than `k` or smaller than `1`. If `NULL`
#' all probabilities of lengths \eqn{s=1,2, \dots, k} are calculated.
#'
#' @returns A list of length `k` with number of possibilities to observe a
#' sequence of misses of length \eqn{s=1, 2, ..., k}. First list element gives
#' the number of possibilities for \eqn{s=1}, second list element for \eqn{2}
#' and so on. Possibilities are of type `bigz` of the `gmp` package. See `?bigz`
#' for further details.
#'
#' @examples
#' freq_seq_misses(23,8)
#'
#' @export

freq_seq_misses <- function(n, k, s=NULL) {

  if(!is.null(s)){

    if(s > k | s < 1){

    stop("s needs to be >0 and <=k.")

    }else{
      (n - k + 1) * gmp::chooseZ(n - s - 1, k - s)
    }

  }else{
    lapply(1:k, function(s) {
      (n - k + 1) * gmp::chooseZ(n - s - 1, k - s)
    })
  }

}
