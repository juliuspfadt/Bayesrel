
# Diagnostics of the fitted factor model underlying the omega coefficients.
# lavaan returns point estimates even when the optimizer did not converge, and when the
# information matrix cannot be inverted the Wald interval of a := parameter collapses onto
# the point estimate (se = 0) instead of being flagged. Both states have to be reported,
# otherwise a coefficient that carries no information reads like an ordinary estimate.
lavCheckFit <- function(fit, sts, labels) {

  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  admissible <- isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))

  se <- sts$se[sts$label %in% labels]
  se.available <- length(se) > 0 && all(is.finite(se)) && all(se > 0)

  return(list(converged = converged, admissible = admissible, se.available = se.available))
}


omegaMultiFreq <- function(data, n.factors, interval, fiml, model, model.type, fit.measures) {

  k <- ncol(data)
  model_opts <- indexMatrix(model, k, n.factors, colnames(data))

  if (model.type == "second-order") {

    if (is.null(model)) {
      modfile <- lavMultiFileSeco(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFileSecoSyntax(k, n.factors, model, colnames(data))

    }

    if (fiml) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = FALSE, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = FALSE)
    }

    sts <- lavaan::parameterestimates(fit, level = interval, standardized = TRUE)
    gloads <- sts$std.all[1:n.factors]
    lmat <- matrix(0, k, n.factors)
    lmat[model_opts$imat] <- sts$std.all[sts$op == "=~" & sts$lhs != "g"]
    theta <- sts$std.all[sts$op == "~~" & sts$lhs %in% colnames(data)]
    psi <- sts$std.all[sts$op == "~~" & sts$lhs %in% modfile$factor_names]

  } else if (model.type == "bi-factor") { # model.type is bifactor
    if (is.null(model)) {

      modfile <- lavMultiFileBif(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified

      if (any(rowSums(model_opts$imat) > 1)) {
        stop("Crossloadings cannot be specified with the bi-factor model.")
      }

      modfile <- lavMultiFileBifSyntax(k, n.factors, model, colnames(data))
    }

    if (fiml) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = TRUE, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = TRUE)
    }

    sts <- lavaan::parameterestimates(fit, level = interval, standardized = TRUE)
    gloads <- sts$std.all[1:k]
    lmat <- matrix(0, k, n.factors)
    lmat[model_opts$imat] <- sts$std.all[(k + 1):(2 * k)]
    psi <- sts$std.all[(2 * k + 2):(2 * k + 1 + n.factors)]
    theta <- sts$std.all[(2 * k + 2 + n.factors):(3 * k + 1 + n.factors)]

  } else if (model.type == "correlated") {

    if (is.null(model)) {
      modfile <- lavMultiFileCorr(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFileCorrSyntax(k, n.factors, model, colnames(data))

    }

    if (fiml) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = FALSE, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = FALSE)
    }

    sts <- lavaan::parameterestimates(fit, level = interval, standardized = TRUE)
    lmat <- matrix(0, k, n.factors)
    lmat[model_opts$imat] <- sts$std.all[sts$op == "=~"]
    theta <- sts$std.all[sts$op == "~~" & sts$lhs %in% colnames(data)]
    psi <- sts$std.all[sts$op == "~~" & sts$lhs %in% modfile$factor_names]

  } else {
    stop("Invalid model type specified.")
  }


  labels <- if (model.type == "correlated") "omega_t" else c("omega_t", "omega_h")
  diagnostics <- lavCheckFit(fit, sts, labels)

  # the fit measures are a separate piece of output: when they cannot be computed, because
  # lavaan refuses them for a model that did not converge, the coefficients themselves are
  # still returned rather than the whole call failing
  if (fit.measures) {
    modfile$fit.measures <- tryCatch(lavaan::fitmeasures(fit), error = function(e) NULL)
    modfile$srmr.summary <- tryCatch(
      if (fiml) lavaan::lavResiduals(fit)$summary["total"] else lavaan::lavResiduals(fit)$summary["cov"],
      error = function(e) NULL)
  }

  if (model.type != "correlated") {
    return(list(omhmean = sts$est[sts$label == "omega_h"], omtmean = sts$est[sts$label == "omega_t"],
                omhlow = sts$ci.lower[sts$label == "omega_h"], omhup = sts$ci.upper[sts$label == "omega_h"],
                omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"],
                lambda = lmat, gloads = gloads, theta = theta, psi = psi,
                modfile = modfile, diagnostics = diagnostics, fit = fit))
  } else {
    return(list(omtmean = sts$est[sts$label == "omega_t"],
                omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"],
                lambda = lmat, theta = theta, psi = psi,
                modfile = modfile, diagnostics = diagnostics, fit = fit))
  }

}
