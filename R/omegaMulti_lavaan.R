
omegaMulti_F <- function(data, n.factors, interval, pairwise) {

  k <- ncol(data)
  modfile <- lavMultiFile_seco(k, n.factors)
  mod <- modfile$model
  colnames(data) <- modfile$names

  if (pairwise) {
    fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F, missing = "ML")
  } else {
    fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F)
  }

  sts <- lavaan::parameterestimates(fit, level = interval)

  return(list(omhmean = sts$est[sts$label == "omega_h"], omtmean = sts$est[sts$label == "omega_t"],
         omhlow = sts$ci.lower[sts$label == "omega_h"], omhup = sts$ci.upper[sts$label == "omega_h"],
         omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"]))
}
