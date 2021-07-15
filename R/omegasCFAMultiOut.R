


omegasCFAMultiOut <- function(data, n.factors, interval, pairwise, model, model.type, fit.measures) {

  out <- list()
  om_out <- omegaMultiF(data, n.factors, interval, pairwise, model, model.type, fit.measures)
  out$omega_t$est <- om_out$omtmean
  out$omega_t$conf <- c(om_out$omtlow, om_out$omtup)

  out$omega_h$est <- om_out$omhmean
  out$omega_h$conf <- c(om_out$omhlow, om_out$omhup)

  out$model <- om_out$modfile
  if (fit.measures) {
    out$fit.measures <- om_out$fit.measures
  }

  return(out)
}