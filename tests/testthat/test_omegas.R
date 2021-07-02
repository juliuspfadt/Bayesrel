

test_that("Bayesian omegas are correct, missing pairwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2)

  expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
               c(0.8635384, 0.8434674, 0.8790391, 0.6422421, 0.5970552, 0.6924216), tolerance = 1e-3)


})


test_that("Frequentist omegas are correct with higher order model, missing pairwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5)

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.864759, 0.8466538, 0.8828642, 0.643837, 0.5943180, 0.6933561), tolerance = 1e-3)


})

test_that("Frequentist omegas are correct with bifactor model, missing listwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5, model.type = "bifactor", missing = "listwise")

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.8718139, 0.8518070, 0.8918207, 0.6309915, 0.5792708, 0.6827122), tolerance = 1e-3)

})
