
#'
#'@export
bomegas <- function(
  data,
  n.factors,
  model = "balanced",
  n.iter = 2e3,
  n.burnin = 200,
  n.chains = 3,
  thin = 1,
  interval = .95,
  missing = "pairwise",
  callback = function(){}
) {

  any_missings <- FALSE
  pairwise <- FALSE
  complete_cases <- nrow(data)
  if (any(is.na(data))) {
    any_missings <- TRUE
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = TRUE)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      complete_cases <- ncomp
    } else { # missing = pairwise
      pairwise <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  sum_res <- bomegas_multi_out(data, n.factors, n.iter, n.burnin, thin, n.chains,
                               interval, model, pairwise, callback)

  sum_res$complete.cases <- complete.cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$pairwise <- pairwise
  class(sum_res) <- 'bomegas'

  return(sum_res)
}
