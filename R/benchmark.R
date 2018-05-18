# library(microbenchmark)
# 
# p <- 25
# a <- matrix(runif(p^2), p, p)
# Sigma <- t(a) %*% a
# xx <- MASS::mvrnorm(n = 1000, rep(0, p), Sigma)
# 
# timin <- microbenchmark(
# 	covSamp(xx, 1e3, 1e2),
# 	covSamp2(xx, 1e3, 1e2),
# 	times = 2e1
# )
# 
# # approx 2 times faster on my machine
# # note that this should become more as the number of iterations increases
# print(timin, unit = "relative") 
# 
# cv <- cov(xx)
# set.seed(1)
# e0 <- LaplacesDemon::rinvwishart(p+2, cv)
# set.seed(1)
# e1 <- rinvwishart(p+2, chol(Matrix::chol2inv(chol(cv))))
# all.equal(e0, e1) # result is equal
# 
# 
# 
