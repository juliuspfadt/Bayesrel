# the basic functions for calculating and bootstrapping the internal consistency estimates

#######  measures functions ##########

applyalpha <- function(M, callback = function(){}){
  a <- alphaArma(M)
  callback()
  return(a)
}

applylambda2 <- function(M, callback = function(){}){
  lambda2 <- l2Arma(M)
  callback()
  return(lambda2)
}

applylambda6 <- function(M, callback = function(){}){

  lambda6 <- try(l6Arma(M), silent = TRUE)
  if (inherits(lambda6, "try-error")) {
    lambda6 <- NaN
    warning("singular bootstrapped covariance matrices encountered when computing lambda6")
  }

  callback()
  return(lambda6)
}

applyomegaPFA <- function(m, callback = function(){}, loadings = FALSE){

  f <- try(pfaArma(m), silent = TRUE)
  if (inherits(f, "try-error")) {
    om <- NaN
    l_fa <- NaN
    warning("singular bootstrapped covariance matrices encountered when computing omega")
  } else {
    l_fa <- f$loadings
    er_fa <- f$err_var
    om <- sum(l_fa)^2 / (sum(l_fa)^2 + sum(er_fa))
    if (om < 0 || om > 1 || is.na(om))
      om <- NaN
  }

  callback()

  if (loadings) {
    return(list(om = om, loadings = l_fa))
  } else {
    return(om)
  }
}

# old functions without Cpp

applyGlbNoCpp <- function(M, callback = function(){}){
  glbObj <- psych::glb.fa(M)
  glb <- glbObj$glb
  callback()
  return(glb)
}

applyalphaNoCpp <- function(M, callback = function(){}){
  p <- ncol(M)
  a <- (p / (p - 1)) * (1 - (sum(diag((M))) / sum(M)))
  callback()
  return(a)
}
applylambda2NoCpp <- function(M, callback = function(){}){
  p <- ncol(M)
  M0 <- M
  diag(M0) <- 0
  lambda2 <- (sum(M0) + sqrt(p / (p - 1) * sum(M0^2))) / sum(M)
  callback()
  return(lambda2)
}

applylambda4NoCpp <- function(M, callback = function(){}){
  if (ncol(M) < 15) {
    out <- MaxSplitExhaustive(M)
  } else {
    l4 <- quant.lambda4(M)
    out <- quantile(l4, prob = 1)
  }
  callback()
  return(out)
}

applylambda6NoCpp <- function(M, callback = function(){}){
  smc <- trySmc(M)
  if (inherits(smc, "try-error") || anyNA(smc)) {
    lambda6 <- NaN
    warning("singular bootstrapped covariance matrices encountered")
  } else {
    lambda6 <- 1 - (sum(1 - (smc)) / sum(cov2cor(M)))
  }
  callback()
  return(lambda6)
}

applyomegaCFAData <- function(data, interval, pairwise, callback = function(){}){
  out <- omegaFreqData(data, interval = interval, omega.int.analytic = TRUE, pairwise = pairwise)
  om <- out$omega
  callback()
  return(om)
}

applyomegaCFACov <- function(cv, interval, omega.int.analytic, pairwise, n.boot){
  data <- MASS::mvrnorm(500, numeric(ncol(cv)), cv)
  out <- omegaFreqData(data, interval, omega.int.analytic, pairwise, n.boot)
  om <- out$omega
  return(om)
}

applyomegaPFANoCpp <- function(m, callback = function(){}){
  f <- princFac(m)
  l_fa <- f$loadings
  er_fa <- f$err_var
  om <- sum(l_fa)^2 / (sum(l_fa)^2 + sum(er_fa))
  if (om < 0 || om > 1 || is.na(om))
    om <- NaN
  callback()
  return(om)
}
