
omegaMulti_F <- function(data, n.factors, interval, pairwise, model.type) {

  k <- ncol(data)
  if (model.type == "higherorder") {
    modfile <- lavMultiFile_seco(k, n.factors)
    mod <- modfile$model
    colnames(data) <- modfile$names

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F)
    }

  } else { # model.type is bifactor
    modfile <- lavMultiFile_bif(k, n.factors)
    mod <- modfile$model
    colnames(data) <- modfile$names

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = T, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = T)
    }
  }


  sts <- lavaan::parameterestimates(fit, level = interval)

  return(list(omhmean = sts$est[sts$label == "omega_h"], omtmean = sts$est[sts$label == "omega_t"],
         omhlow = sts$ci.lower[sts$label == "omega_h"], omhup = sts$ci.upper[sts$label == "omega_h"],
         omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"]))
}
