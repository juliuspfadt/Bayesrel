# forces a quadratic matrix to be symmetrical
make_symmetric <- function(a, lower.tri = TRUE) {
  if (lower.tri){
    ind <- upper.tri(a)
    a[ind] <- t(a)[ind]
  } else {
    ind <- lower.tri(a)
    a[ind] <- t(a)[ind]
  }
  return(a)
}

# computes alpha analytical interval with given bounds
# source:
# Bonett, D. G., & Wright, T. A. (2015). Cronbach’s alpha reliability:
# Interval estimation, hypothesis testing, and sample size planning. Journal of Organizational Behavior,36(1), 3–15.
#
ciAlpha <- function(palpha, n, V) {
  p <- ncol(V)
  z <- qnorm(1 - palpha/2)
  b <- log(n/ (n - 1))
  j <- matrix(rep(1, p), nrow = p, ncol = 1)
  a0 <- t(j) %*% V %*% j
  t1 <- sum(diag(V))
  t2 <- sum(diag(V %*% V))
  a1 <- a0 ^ 3
  a2 <- a0 * (t2 + t1^2)
  a3 <- 2 * t1 * t(j) %*% (V %*% V) %*% j
  a4 <- 2 * p^2 / (a1 * (p - 1)^2)
  r <- (p/ (p-1)) * (1 - t1 / a0)
  var <- a4 * (a2 - a3) / (n - 3)
  ll <- 1 - exp(log(1 - r) - b + z * sqrt(var / (1 - r)^2))
  ul <- 1 - exp(log(1 - r) - b - z * sqrt(var / (1 - r)^2))
  out <- c(ll, ul)
  return(out)
}

# does quantile averaging and returns 2000 datapoints
quantiles <- function(samp, length_out = 2e3){
  q <- quantile(samp, probs = seq(0, 1, length.out = length_out))
  return(q)
}


se <- function(x) {
  b <- length(x)
  se <- sqrt(1 / (b-1) * sum((x - mean(x))^2))
  se
}


# create lavaan cfa one factor model file from data

lavOneFile <- function(data){
  p <- ncol(data)
  v <- 0
  for (i in 1:p) {
    v[i] <- paste0("x", i)
  }
  v <- paste0(v, collapse = "+")
  mod <- paste0("g=~", v) # dynamic lavaan model file
  mod <- paste0(mod, "; g ~~ 1*g") # fix the factor variance to 1

  # column names specify
  names <- 0
  for (i in 1:p) {
    names[i] <- paste0("x", i)
  }
  return(list(names = names, model = mod))
}


# calculate omega from loadings and residual (error variances)

omegaBasic <- function(l, e) {
  o <- sum(l)^2 / (sum(l)^2 + sum(e))
  return(o)
}

# calculate the kolomogorov smirnov distances between some samples and the original sample
KSTestStatistic <- function(x, y) {
  t <- stats::ks.test(x, y)
  return(t$statistic)
}

# calculate the kublack leibler distance between two samples
KLDStatistic <- function(x, y) {
  # transform the samples to PDFs:
  xdf <- get_approx_density(x)
  ydf <- get_approx_density(y)

  xx <- seq(0, 1, length.out = 1e3)
  t <- LaplacesDemon::KLD(xdf(xx), ydf(xx))
  return(t$sum.KLD.py.px)
}

hpdHelp <- function(x) {
  x <- coda::as.mcmc(x)
  return(coda::HPDinterval(x))
}


# create covariance matrix
createUnidimCovMat <- function(avg, p) {
  mean_cor <- 1
  counter <- 1
  while (mean_cor < (avg - .001) || mean_cor > (avg + .001)) {
    mlam <- avg * 3 + .02
    vlam <- avg * 2
    lam_true <- abs(rnorm(p, mlam, vlam))
    psi_true <- 1/rgamma(p, 10, 10)
    loading <- matrix(lam_true, nrow = p)
    psi_m <- diag(1, nrow = p)
    diag(psi_m) <- psi_true
    tmpCov <- make_symmetric(loading %*% 1 %*% t(loading) + psi_m)
    cormat <- cov2cor(tmpCov)
    mean_cor <- (sum(cormat) - p) / (p*p - p)
    counter <- counter + 1
    if (counter == 1e4)
      return(print("solution has not been found"))
  }
  return(tmpCov)
}

try_smc <- function(M) {
  return(try(1 - 1 / diag(solve(cov2cor(M))), silent = TRUE))
}

try_psd <- function(M) {
  R <- cov2cor(M)
  return(try(eigen(R, only.values = TRUE)$values, silent = TRUE))
}


get_approx_density <- function(x) {
  d <- density(x, n = 2^12)
  f <- approxfun(d$x, d$y, yleft = 0, yright = 0)
  c <- integrate(f, 0, 1)$value
  return(
    function(x) {
      return(f(x) / c)
    }
  )
}

omegas_seco <- function(lambda, beta, theta, psi) {
  gl <- lambda %*% beta
  sl <- lambda %*% sqrt(psi)
  omh <- sum(gl %*% t(gl)) / (sum(gl %*% t(gl)) + sum(sl %*% t(sl)) + sum(theta))
  omt <- (sum(gl %*% t(gl)) + sum(sl %*% t(sl))) / (sum(gl %*% t(gl)) + sum(sl %*% t(sl)) + sum(theta))

  return(c(
    omh, omt
  ))
}

# multivariate normal data with matrix of means and V
genNormData_tweak <- function(n, m, Sigma){
  p <- ncol(Sigma)
  randomData <- matrix(rnorm(n*p), n, p)
  cc <- chol(Sigma)
  out <- randomData %*% cc
  out <- out + matrix(m, nrow = n, ncol = ncol(Sigma), byrow = FALSE)
  return(out)
}

# get model implied covariance matrix
implCov_multi <- function(s, b, theta, psi) {
  i <- diag(nrow(b))
  ib_inv <- solve(i-b)
  out <- s %*% ib_inv %*% psi %*% t(ib_inv) %*% t(s) + theta
  out <- Bayesrel:::make_symmetric(out)
  return(out)
}

# generate balanced and simple model file for lavaan
lavMultiFile_seco <- function(k, ns) {

  beta_names <- paste("gl", 1:ns, sep = "")
  lambda_names <- matrix(paste("sl", 1:k, sep = ""), k/ns, ns)
  names <- matrix(paste("x", 1:k, sep = ""), k/ns, ns)
  theta_names <- paste("e", 1:k, sep = "")
  psi_names <- paste("ps", 1:ns, sep = "")

  gen <- paste(paste(beta_names, "*", "s", 1:ns, sep = ""), collapse = " + " )
  gen <- paste0("g =~ ", gen)

  mod <- NULL
  for (i in 1:ns) {
    mod <- paste0(mod, "s", i, " =~ ", paste(paste(lambda_names[, i], "*", names[, i], sep = ""),
                                             collapse = " + "), "\n")
  }
  mod <- paste0(gen, " \n", mod)

  # extra stuff for omega calc
  mod <- paste0(mod, "g ~~ 1*g\n")
  mod <- paste0(mod, paste(paste(c(names), " ~~ ", theta_names, "*",
                                 c(names), sep = ""), collapse = "\n"), "\n")

  # define extra parameters
  gl <- NULL
  sl <- NULL
  for (i in 1:ns) {
    gl <- c(gl, paste0(beta_names[i] , "*", lambda_names[, i]))
    sl <- paste0(sl, "ssum", i, " := ", paste0(lambda_names[, i], collapse = " + "), "\n")
  }
  gl <- paste0(gl, collapse = "+")
  sum_gl <- paste("g_loading :=", gl, "\n")


  sum_sl <- paste0("spec_loading := ", paste(paste0("ssum", 1:ns, "^2"), collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_h <- "omega_h := (g_loading^2) / (g_loading^2 + spec_loading + residual_var) \n"
  omega_t <- "omega_t := (g_loading^2 + spec_loading) / (g_loading^2 + spec_loading + residual_var) \n"
  mod <- paste0(mod, sl, sum_gl, sum_sl, sum_errs, omega_h, omega_t)

  return(out <- list(names = names, model = mod))
}

