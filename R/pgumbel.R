#' Cumulative distribution function of the Gumbel distribution
#'
#' Cumulative distribution function of the Gumbel distribution
#' \deqn{\exp^{-\exp^\frac{\mu - k}{s}}.}
#'
#' @param k Realisation of the support.
#' @param mu The location parameter.
#' @param s The scale parameter.
#' @returns The p-values for each gene to observe its biggest gap.
#' @examples
#' pgumbel(10,5,3)
#' @export

pgumbel <- function(k, mu, s){
  exp(
    - exp(
      (mu - k) / s
    )
  )
}

