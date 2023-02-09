
tol <- ifelse(testthat:::on_cran(), 1e-2, 1e-3)

test_that("Bayesian omegas are correct, missing pairwise, and fit indices are good", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2, model.type = "second-order")

  expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
               c(0.8633604, 0.8448201, 0.8790073, 0.6399786, 0.5956023, 0.6922816),
               tolerance = tol)

  ff <- secoFit(ee, upps, ppc = FALSE, cutoff = .06)
  expect_equal(as.numeric(c(ff$LR, ff$srmr_pointEst, ff$rmsea_pointEst, ff$rmsea_ci, ff$p_rmsea)),
               c(409.03892661, 0.07028219, 0.05616901, 0.05424941, 0.05760252,
                 1.00000000), tolerance = tol)

})

test_that("Bayesian omegas are correct, missing listwise, model sytnax specified", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  mod <- "
  f1 =~ U17_r+U22_r+U29_r+U34_r
  f2 =~ U4+U14+U19+U27
  f3 =~ U6 +U16+U28+U48
  f4 =~ U23_r +U31_r +U36_r +U46_r
  f5 =~ U10_r +U20_r +U35_r +U52_r
  "
  ee <- Bayesrel::bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2,
                          missing = "listwise", model = mod, model.type = "second-order")

  expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
               c(0.8624262, 0.8451732, 0.8801895, 0.6402192, 0.5909983, 0.6860989), tolerance = tol)

})


test_that("Frequentist omegas are correct with higher order model, missing pairwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5)

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.864759, 0.8466538, 0.8828642, 0.643837, 0.5943180, 0.6933561), tolerance = tol)


})

test_that("Frequentist omegas are correct with bifactor model, missing listwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5, model.type = "bi-factor", missing = "listwise")

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.8718139, 0.8518070, 0.8918207, 0.6309915, 0.5792708, 0.6827122), tolerance = tol)

})


test_that("Frequentist omegas are correct with bifactor model, missing listwise, model sytnax specified,
          first and last item switched", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  mod <- "
  f1 =~ U52_r+U22_r+U29_r+U34_r
  f2 =~ U4+U14+U19+U27
  f3 =~ U6 +U16+U28+U48
  f4 =~ U23_r +U31_r +U36_r +U46_r
  f5 =~ U10_r +U20_r +U35_r + U17_r
  "
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5, model.type = "bi-factor", missing = "listwise",
                            model = mod)

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.8702046, 0.8499908, 0.8904183, 0.6275952, 0.5760795, 0.6791109), tolerance = tol)

})


test_that("Bayesian omegas are correct with altered prior hyperparameters", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2,
                          a0 = 6, b0 = 10, l0 = 1, c0 = 10, d0 = 6, beta0 = 2, p0 = 12, R0 = 5,
                          model.type = "second-order")

  expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
               c(0.8476525, 0.8275225, 0.8644323, 0.6255007, 0.5799303, 0.6704586),
               tolerance = tol)

})

test_that("Bayesian omegas are correct with cross loadings", {

  data(upps, package = "Bayesrel")
  mod <- "
  f1 =~ U17_r + U22_r + U29_r + U34_r + U19
  f2 =~ U4 + U14 + U19 + U27
  f3 =~ U6 + U16 + U28 + U48 + U29_r
  f4 =~ U23_r + U31_r + U36_r + U46_r + U34_r
  f5 =~ U10_r + U20_r + U35_r + U52_r + U19
  "
  set.seed(1234)
  ee <- Bayesrel::bomegas(upps, n.iter = 200, n.burnin = 50, n.chains = 2,
                          model = mod, missing = "listwise", param.out = TRUE, model.type = "second-order")
  ll <- apply(ee$model$lambda, c(3, 4), mean)

  expect_equal(ll,
               matrix(c(-0.5236519, -0.5287243, -0.4082028, -0.628011, 0, 0, 0.07909764, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.003327511, 0.004258746, -0.0007414194,
                        0.001434998, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001336209, 0, 0, 0, 0,
                        0, -0.0008860797, 0.003827371, 0.002945267, -0.002850434, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0.005032167, 0, 0, 0, 0, 0, 0, 0, 0, -0.004571036, 0.001070471, -0.002486544,
                        0.003125839, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.002564826, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0.003184005, 0.003274632, 0.003072611, -0.004328613), 20, 5),
               tolerance = tol)

})


test_that("Frequentist omegas are correct with cross loadings", {

  data(upps, package = "Bayesrel")
  mod <- "
  f1 =~ U17_r + U22_r + U29_r + U34_r + U19
  f2 =~ U4 + U14 + U19 + U27
  f3 =~ U6 + U16 + U28 + U48 + U29_r
  f4 =~ U23_r + U31_r + U36_r + U46_r + U34_r
  f5 =~ U10_r + U20_r + U35_r + U52_r + U19
  "
  ee <- Bayesrel::omegasCFA(upps, model = mod, missing = "listwise")
  expect_equal(as.numeric(c(unlist(ee[[1]]), unlist(ee[[2]]))),
               c(0.8648853, 0.8466552, 0.8831155, 0.6464642, 0.5962111, 0.6967174),
               tolerance = tol)

})
