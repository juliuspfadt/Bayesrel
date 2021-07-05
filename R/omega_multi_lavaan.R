
omegaMulti_F <- function(data, n.factors, interval, pairwise, model, model.type) {

  k <- ncol(data)
  if (model.type == "higher.order") {

    if (model == "balanced") {
      modfile <- lavMultiFile_seco(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFile_seco_syntax(k, n.factors, model, colnames(data))

      if (!modfile$colnames) { # name the variables in the dataset:
        names <- unlist(modfile$names)
        inds <- as.numeric(unlist(regmatches(names, gregexpr("[[:digit:]]+", names))))
        colnames(data)[inds] <- names
      }
    }

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F)
    }


  } else if (model.type == "bi.factor") { # model.type is bifactor
    if (model == "balanced") {

      modfile <- lavMultiFile_bif(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFile_bif_syntax(k, n.factors, model, colnames(data))

      if (!modfile$colnames) { # name the variables in the dataset:
        names <- unlist(modfile$names)
        inds <- as.numeric(unlist(regmatches(names, gregexpr("[[:digit:]]+", names))))
        colnames(data)[inds] <- names
      }
    }

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = T, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = T)
    }
  }


  sts <- lavaan::parameterestimates(fit, level = interval)

  return(list(omhmean = sts$est[sts$label == "omega_h"], omtmean = sts$est[sts$label == "omega_t"],
         omhlow = sts$ci.lower[sts$label == "omega_h"], omhup = sts$ci.upper[sts$label == "omega_h"],
         omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"],
         modfile = modfile))
}
