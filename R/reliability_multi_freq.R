


#'
#'@export
omegasCFA <- function(
  data,
  n.factors,
  model = "balanced",
  model.type = "higherorder",
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
      complete.cases <- ncomp
    } else { # missing = pairwise
      pairwise <- TRUE
      complete.cases <- nrow(data)
    }
  }

  data <- scale(data, scale = FALSE)

  sum_res <- omegasCFA_multi_out(data, n.factors, interval, pairwise, model.type)

  sum_res$complete.cases <- complete.cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$pairwise <- pairwise
  class(sum_res) <- 'omegasCFA'

  return(sum_res)
}
