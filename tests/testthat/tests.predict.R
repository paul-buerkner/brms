test_that("predict for location shift models runs without errors", {
  ns <- 30
  nobs <- 10
  draws <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                sigma = rchisq(ns, 3), nu = rgamma(ns, 4),
                nsamples = ns)
  i <- sample(nobs, 1)
  
  draws$f$link <- "identity"
  pred <- predict_gaussian(i, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$f$link <- "log"
  pred <- predict_student(i, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$f$link <- "inverse"
  pred <- predict_cauchy(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for lognormal and exgaussian models runs without errors", {
  ns <- 50
  nobs <- 2
  draws <- list(sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
                mu = matrix(rnorm(ns * nobs), ncol = nobs),
                f = lognormal(), nsamples = ns)
  pred <- predict_lognormal(1, draws = draws)
  expect_equal(length(pred), ns)
  pred <- predict_exgaussian(1, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for aysm_laplace models runs without errors", {
  ns <- 50
  draws <- list(sigma = rchisq(ns, 3), quantile = rbeta(ns, 2, 1),
                mu = matrix(rnorm(ns*2), ncol = 2),
                f = asym_laplace(), nsamples = ns)
  pred <- brms:::predict_asym_laplace(1, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), 
                dim = c(3, 3, 10))
  draws <- list(mu = matrix(rnorm(ns*nobs), ncol = nobs),
            Sigma = aperm(Sigma, c(3, 1, 2)), 
            nu = matrix(rgamma(ns, 5)),
            nsamples = ns)
  draws$data <- list(N = nobs, N_trait = ncols)
  draws$f$link <- "identity"
  
  pred <- predict_gaussian_mv(1, draws = draws)
  expect_equal(dim(pred), c(ns, nvars))
  
  pred <- predict_student_mv(2, draws = draws)
  expect_equal(dim(pred), c(ns, nvars))
  
  pred <- predict_cauchy_mv(3, draws = draws)
  expect_equal(dim(pred), c(ns, nvars))
})

test_that("predict for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  draws <- list(mu = matrix(rnorm(ns*nobs), ncol = nobs),
                sigma = matrix(rchisq(ns, 3)),
                nu = matrix(rgamma(ns, 5)),
                ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
                ma = matrix(rnorm(ns, 0.2, 1), ncol = 1),
                nsamples = ns)
  draws$data <- list(begin_tg = c(1, 5, 12), nobs_tg = c(4, 7, 3),
               se2 = rgamma(ns, 10))
  
  draws$f$link <- "inverse"
  pred <- predict_gaussian_cov(1, draws = draws)
  expect_equal(length(pred), ns * 4)
  
  draws$f$link <- "log"
  pred <- predict_student_cov(2, draws = draws[-4])
  expect_equal(length(pred), ns * 7)
  
  draws$f$link <- "identity" 
  pred <- predict_cauchy_cov(3, draws = draws[-5])
  expect_equal(length(pred), ns * 3)
})

test_that("predict for 'cor_fixed' models runs without errors", {
  draws <- list(mu = matrix(rnorm(30), nrow = 3),
                nu = matrix(rep(2, 3)), nsamples = 3)
  draws$data <- list(V = diag(10))
  draws$f$link <- "identity"
  pred <- predict_gaussian_fixed(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
  pred <- predict_student_fixed(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
  pred <- predict_cauchy_fixed(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
})

test_that("predict for count and survival models runs without errors", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  draws <- list(mu = matrix(rnorm(ns*nobs), ncol = nobs),
                shape = rgamma(ns, 4), xi = 0,
                nsamples = ns)
  draws$nu <- draws$sigma <- draws$shape + 1
  draws$data <- list(trials = trials)
  i <- sample(nobs, 1)
  
  draws$f$link <- "cloglog"
  pred <- predict_binomial(i, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$f$link <- "log"
  pred <- predict_poisson(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_negbinomial(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_geometric(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_exponential(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_gamma(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_weibull(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_frechet(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_inverse.gaussian(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_gen_extreme_value(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for bernoulli and beta models works correctly", {
  ns <- 17
  nobs <- 10
  draws <- list(mu = matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs),
                phi = rgamma(ns, 4), nsamples = ns)
  i <- sample(1:nobs, 1)
  draws$f$link <- "logit"
  
  pred <- predict_bernoulli(i, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$data <- list(N_trait = nobs)
  pred <- predict_bernoulli(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_beta(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  draws <- list(mu = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
                kappa = matrix(rgamma(ns, 4)), nsamples = ns)
  draws$f$link <- "tan_half"
  i <- sample(seq_len(nobs), 1)
  pred <- predict_von_mises(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  draws <- list(mu = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
                shape = rgamma(ns, 4), phi = rgamma(ns, 1),
                zi = rbeta(ns, 1, 1), coi = rbeta(ns, 5, 7), 
                nsamples = ns)
  draws$hu <- draws$zoi <- draws$zi
  draws$data <- list(N_trait = nobs, trials = trials)
  draws$f$link <- "log"
  
  pred <- predict_hurdle_poisson(1, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_hurdle_negbinomial(2, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_hurdle_gamma(5, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_poisson(3, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_negbinomial(6, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$f$link <- "logit"
  pred <- predict_zero_inflated_binomial(4, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_beta(8, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_one_inflated_beta(7, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for categorical and ordinal models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 4
  draws <- list(mu = array(rnorm(ns*nobs), dim = c(ns, nobs, ncat)),
                nsamples = ns)
  draws$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  
  draws$f$link <- "logit"
  pred <- sapply(1:nobs, predict_categorical, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_cumulative, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_sratio, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_cratio, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_acat, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$f$link <- "probit"
  pred <- sapply(1:nobs, predict_acat, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("truncated predict run without errors", {
  ns <- 30
  nobs <- 15
  draws <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                sigma = rchisq(ns, 3), nsamples = ns)
  
  draws$f$link <- "identity"
  draws$data <- list(lb = sample(-(4:7), nobs, TRUE))
  pred <- sapply(1:nobs, predict_gaussian, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$f$link <- "log"
  draws$data <- list(ub = sample(70:80, nobs, TRUE))
  pred <- sapply(1:nobs, predict_poisson, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$data <- list(lb = rep(0, nobs), ub = sample(70:75, nobs, TRUE))
  pred <- sapply(1:nobs, predict_poisson, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("predict for the wiener diffusion model runs without errors", {
  ns <- 5
  nobs <- 3
  draws <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                bs = matrix(rchisq(ns, 3)), ndt = matrix(rep(0.5, ns)),
                bias = matrix(rbeta(ns, 1, 1)))
  draws$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
  draws$f$link <- "identity"
  i <- sample(1:nobs, 1)
  expect_equal(nrow(predict_wiener(i, draws)), ns)
})
