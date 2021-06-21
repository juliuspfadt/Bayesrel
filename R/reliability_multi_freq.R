


#'
#'@export
omegasCFA <- function(
  data,
  n.factors,
  model = "balanced",
  interval = .95,
  missing = "pairwise") {

  any_missings <- FALSE
  pairwise <- FALSE
  if (any(is.na(data))) {
    any_missings <- TRUE
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = TRUE)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      sum_res$complete.cases <- ncomp
    } else { # missing = pairwise
      pairwise <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  sum_res <- omegasCFA_multi_out(data, n.factors, interval, pairwise)


  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$pairwise <- pairwise
  class(sum_res) <- 'omegasCFA'

  return(sum_res)
}
