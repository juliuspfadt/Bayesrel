



omegasCFA_multi_out <- function(data, n.factors, interval, pairwise) {

  out <- list()
  om_out <- omegaMulti_F(data, n.factors, interval, pairwise)
  out$omega_t$est <- om_out$omtmean
  out$omega_t$conf <- c(om_out$omtlow, om_out$omtup)

  out$omega_h$est <- om_out$omhmean
  out$omega_h$conf <- c(om_out$omhlow, om_out$omhup)

  return(out)
}
