


omegasCFAMultiOut <- function(data, n.factors, interval, fiml, model, model.type, fit.measures) {

  out <- list()
  om_out <- omegaMultiFreq(data, n.factors, interval, fiml, model, model.type, fit.measures)

  # without usable standard errors the Wald interval collapses onto the point estimate;
  # reporting NA is honest, a zero-width interval would claim perfect precision
  conf <- function(low, up) if (om_out$diagnostics$se.available) c(low, up) else c(NA_real_, NA_real_)

  out$omega_t$est <- om_out$omtmean
  out$omega_t$conf <- conf(om_out$omtlow, om_out$omtup)

  out$omega_h$est <- om_out$omhmean
  out$omega_h$conf <- conf(om_out$omhlow, om_out$omhup)

  out$loadings$specific <- om_out$lambda
  out$residuals$specific <- om_out$theta
  out$loadings$general <- om_out$gloads
  out$residuals$general <- om_out$psi

  out$model <- om_out$modfile
  if (fit.measures) {
    out$fit.measures <- om_out$fit.measures
  }

  out$diagnostics <- om_out$diagnostics
  out$lavaan.fit <- om_out$fit

  return(out)
}
