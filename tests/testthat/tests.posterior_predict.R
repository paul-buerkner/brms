context("Tests for posterior_predict helper functions")

test_that("posterior_predict for location shift models runs without errors", {
  ns <- 30
  nobs <- 10
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3), nu = rgamma(ns, 4)
  )
  i <- sample(nobs, 1)

  pred <- brms:::posterior_predict_gaussian(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_student(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for various skewed models runs without errors", {
  ns <- 50
  nobs <- 2
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    alpha = rnorm(ns), ndt = 1
  )
  pred <- brms:::posterior_predict_lognormal(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_shifted_lognormal(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_exgaussian(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_skew_normal(1, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for aysm_laplace models runs without errors", {
  ns <- 50
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    sigma = rchisq(ns, 3),
    quantile = rbeta(ns, 2, 1),
    mu = matrix(rnorm(ns*2), ncol = 2),
    zi = rbeta(ns, 10, 10)
  )
  pred <- brms:::posterior_predict_asym_laplace(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_zero_inflated_asym_laplace(1, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), dim = c(3, 3, 10))
  prep <- structure(list(ndraws = ns), class = "mvbrmsprep")
  prep$mvpars <- list(
    Mu = array(rnorm(ns*nobs*nvars), dim = c(ns, nobs, nvars)),
    Sigma = aperm(Sigma, c(3, 1, 2))
  )
  prep$dpars <- list(nu = rgamma(ns, 5))
  prep$data <- list(N = nobs, N_trait = ncols)

  pred <- brms:::posterior_predict_gaussian_mv(1, prep = prep)
  expect_equal(dim(pred), c(ns, nvars))

  pred <- brms:::posterior_predict_student_mv(2, prep = prep)
  expect_equal(dim(pred), c(ns, nvars))
})

test_that("posterior_predict for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns*nobs), ncol = nobs),
    sigma = rchisq(ns, 3),
    nu = rgamma(ns, 5)
  )
  prep$ac <- list(
    ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
    ma = matrix(rnorm(ns, 0.2, 1), ncol = 1),
    begin_tg = c(1, 5, 12), end_tg = c(4, 11, 15)
  )
  prep$data <- list(se = rgamma(ns, 10))

  prep$family$fun <- "gaussian_time"
  pred <- brms:::posterior_predict_gaussian_time(1, prep = prep)
  expect_equal(length(pred), ns * 4)

  prep$family$fun <- "student_time"
  pred <- brms:::posterior_predict_student_time(2, prep = prep)
  expect_equal(length(pred), ns * 7)
})

test_that("loglik for SAR models runs without errors", {
  ns = 3
  prep <- structure(list(ndraws = ns, nobs = 10), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(30), nrow = ns),
    nu = rep(2, ns),
    sigma = rep(10, ns)
  )
  prep$ac <- list(lagsar = matrix(c(0.3, 0.5, 0.7)), Msar = diag(10))

  pred <- brms:::posterior_predict_gaussian_lagsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::posterior_predict_student_lagsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))

  prep$ac$errorsar <- prep$ac$lagsar
  prep$ac$lagsar <- NULL
  pred <- brms:::posterior_predict_gaussian_errorsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::posterior_predict_student_errorsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))
})

test_that("posterior_predict for FCOR models runs without errors", {
  ns <- 3
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(nobs * ns), nrow = ns),
    sigma = rep(1, ns), nu = rep(2, ns)
  )
  prep$ac <- list(Mfcor = diag(nobs))
  pred <- brms:::posterior_predict_gaussian_fcor(1, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
  pred <- brms:::posterior_predict_student_fcor(1, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_predict for count and survival models runs without errors", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns * nobs), ncol = nobs),
    shape = rgamma(ns, 4), xi = 0, phi = rgamma(ns, 1)
  )
  prep$dpars$nu <- prep$dpars$sigma <- prep$dpars$shape + 1
  prep$data <- list(trials = trials)
  i <- sample(nobs, 1)

  prep$dpars$mu <- brms:::inv_cloglog(prep$dpars$eta)
  pred <- brms:::posterior_predict_binomial(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_beta_binomial(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_discrete_weibull(i, prep = prep)
  expect_equal(length(pred), ns)

  prep$dpars$mu <- exp(prep$dpars$eta)
  pred <- brms:::posterior_predict_poisson(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_negbinomial(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_negbinomial2(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_geometric(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_com_poisson(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_exponential(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_gamma(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_frechet(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_inverse.gaussian(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_gen_extreme_value(i, prep = prep)
  expect_equal(length(pred), ns)

  prep$family$link <- "log"
  pred <- brms:::posterior_predict_weibull(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for bernoulli and beta models works correctly", {
  ns <- 17
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs)),
    phi = rgamma(ns, 4)
  )
  i <- sample(1:nobs, 1)

  pred <- brms:::posterior_predict_bernoulli(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_beta(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for xbeta models works correctly", {
  skip_if_not_installed("betareg")

  ns <- 17
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs)),
    phi = rgamma(ns, 4),
    kappa = rexp(ns, 5)
  )
  i <- sample(1:nobs, 1)

  pred <- brms:::posterior_predict_xbeta(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_xbeta(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = 2 * atan(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    kappa = rgamma(ns, 4)
  )
  i <- sample(seq_len(nobs), 1)
  pred <- brms:::posterior_predict_von_mises(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
    shape = rgamma(ns, 4), phi = rgamma(ns, 1),
    zi = rbeta(ns, 1, 1), coi = rbeta(ns, 5, 7)
  )
  prep$dpars$hu <- prep$dpars$zoi <- prep$dpars$zi
  prep$data <- list(trials = trials)

  prep$dpars$mu <- exp(prep$dpars$eta)
  pred <- brms:::posterior_predict_hurdle_poisson(1, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_hurdle_negbinomial(2, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_hurdle_gamma(5, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_poisson(3, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_negbinomial(6, prep = prep)
  expect_equal(length(pred), ns)

  prep$dpars$mu <- brms:::inv_logit(prep$dpars$eta)
  pred <- brms:::posterior_predict_zero_inflated_binomial(4, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_beta_binomial(6, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_beta(8, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_one_inflated_beta(7, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for ordinal models runs without errors", {
  ns <- 50
  nobs <- 8
  nthres <- 3
  ncat <- nthres + 1
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = array(rnorm(ns * nobs), dim = c(ns, nobs)),
    disc = rexp(ns),
    hu = rbeta(ns, 1, 1)
  )
  prep$thres$thres <- array(0, dim = c(ns, nthres))
  prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  prep$family$link <- "logit"

  prep$family$family <- "cumulative"
  pred <- sapply(1:nobs, brms:::posterior_predict_cumulative, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "sratio"
  pred <- sapply(1:nobs, brms:::posterior_predict_sratio, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "cratio"
  pred <- sapply(1:nobs, brms:::posterior_predict_cratio, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "acat"
  pred <- sapply(1:nobs, brms:::posterior_predict_acat, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$link <- "probit"
  pred <- sapply(1:nobs, brms:::posterior_predict_acat, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "hurdle_cumulative"
  pred <- sapply(1:nobs, brms:::posterior_predict_hurdle_cumulative, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_predict for categorical and related models runs without erros", {
  set.seed(1234)
  ns <- 50
  nobs <- 8
  ncat <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu1 = array(rnorm(ns*nobs, 0, 0.1), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs, 0, 0.1), dim = c(ns, nobs))
  )
  prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  prep$family <- categorical()
  prep$refcat <- 1
  pred <- sapply(1:nobs, brms:::posterior_predict_categorical, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$data$trials <- sample(1:20, nobs)
  prep$family <- multinomial()
  pred <- brms:::posterior_predict_multinomial(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))

  prep$data$trials <- sample(1:20, nobs)
  prep$dpars$phi <- rexp(ns, 1)
  prep$family <- dirichlet_multinomial()
  pred <- brms:::posterior_predict_dirichlet_multinomial(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))

  prep$dpars$phi <- rexp(ns, 1)
  prep$family <- dirichlet()
  pred <- brms:::posterior_predict_dirichlet(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))

  prep$family <- brmsfamily("dirichlet2")
  prep$dpars$mu1 <- rexp(ns, 10)
  prep$dpars$mu2 <- rexp(ns, 10)
  prep$dpars$mu3 <- rexp(ns, 10)
  pred <- brms:::posterior_predict_dirichlet2(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))

  prep$family <- brmsfamily("logistic_normal")
  prep$dpars <- list(
    mu2 = rnorm(ns),
    mu3 = rnorm(ns),
    sigma2 = rexp(ns, 10),
    sigma3 = rexp(ns, 10)
  )
  prep$lncor <- rbeta(ns, 2, 1)
  pred <- brms:::posterior_predict_logistic_normal(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))
})

test_that("truncated posterior_predict run without errors", {
  ns <- 30
  nobs <- 15
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3)
  )
  prep$refcat <- 1

  prep$data <- list(lb = sample(-(4:7), nobs, TRUE))
  pred <- sapply(1:nobs, brms:::posterior_predict_gaussian, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$dpars$mu <- exp(prep$dpars$mu)
  prep$data <- list(ub = sample(70:80, nobs, TRUE))
  pred <- sapply(1:nobs, brms:::posterior_predict_poisson, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$data <- list(lb = rep(0, nobs), ub = sample(70:75, nobs, TRUE))
  pred <- sapply(1:nobs, brms:::posterior_predict_poisson, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_predict for the wiener diffusion model runs without errors", {
  skip("skip as long as RWiener fails on R-devel for 3.6.0")
  ns <- 5
  nobs <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    bs = rchisq(ns, 3), ndt = rep(0.5, ns),
    bias = rbeta(ns, 1, 1)
  )
  prep$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
  i <- sample(1:nobs, 1)
  expect_equal(nrow(brms:::posterior_predict_wiener(i, prep)), ns)
})

test_that("posterior_predict_custom runs without errors", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rbeta(ns * nobs * 2, 1, 1), ncol = nobs * 2)
  )
  prep$data <- list(trials = rep(1, nobs))
  prep$family <- custom_family(
    "beta_binomial2", dpars = c("mu", "tau"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "trials[n]"
  )
  posterior_predict_beta_binomial2 <- function(i, prep) {
    mu <- prep$dpars$mu[, i]
    rbinom(prep$ndraws, size = prep$data$trials[i], prob = mu)
  }
  expect_equal(length(brms:::posterior_predict_custom(sample(1:nobs, 1), prep)), ns)
})

make_prep_outcome <- function(
    ns = 120, nobs = 8, seed = 1001, dpars = list(), data = list()) {
  set.seed(seed)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- dpars
  prep$data <- modifyList(
    list(Y = rnorm(nobs), lb = rep(NULL, nobs), ub = rep(NULL, nobs)),
    data
  )
  prep
}

make_prep_gaussian_outcome <- function(ns = 120, nobs = 8) {
  make_prep_outcome(
    ns = ns, nobs = nobs, seed = 1001,
    dpars = list(
      mu = matrix(rnorm(ns * nobs), ncol = nobs),
      sigma = rgamma(ns, shape = 4, rate = 3)
    )
  )
}

make_prep_student_outcome <- function(ns = 120, nobs = 8) {
  make_prep_outcome(
    ns = ns, nobs = nobs, seed = 1002,
    dpars = list(
      mu = matrix(rnorm(ns * nobs), ncol = nobs),
      sigma = rgamma(ns, shape = 4, rate = 3),
      nu = rgamma(ns, shape = 6, rate = 1) + 2
    )
  )
}

make_prep_positive_outcome <- function(ns = 120, nobs = 8) {
  make_prep_outcome(
    ns = ns, nobs = nobs, seed = 1003,
    dpars = list(
      mu = matrix(exp(rnorm(ns * nobs, mean = 0.2, sd = 0.4)), ncol = nobs),
      sigma = rgamma(ns, shape = 4, rate = 3),
      beta = rgamma(ns, shape = 4, rate = 3),
      shape = rgamma(ns, shape = 6, rate = 2) + 0.5,
      ndt = runif(ns, min = 0, max = 0.5),
      alpha = rnorm(ns),
      nu = rgamma(ns, shape = 4, rate = 1) + 2,
      xi = rnorm(ns, sd = 0.3),
      quantile = runif(ns, min = 0.2, max = 0.8),
      phi = rgamma(ns, shape = 5, rate = 1),
      kappa = rgamma(ns, shape = 2, rate = 1),
      zi = rbeta(ns, 1.5, 5)
    ),
    data = list(Y = rgamma(nobs, shape = 2, rate = 1))
  )
}

make_prep_beta_outcome <- function(ns = 140, nobs = 9) {
  make_prep_outcome(
    ns = ns, nobs = nobs, seed = 1004,
    dpars = list(
      mu = matrix(plogis(rnorm(ns * nobs)), ncol = nobs),
      phi = rgamma(ns, shape = 5, rate = 1),
      kappa = rgamma(ns, shape = 2, rate = 1),
      zi = rbeta(ns, 1.5, 5),
      zoi = rbeta(ns, 1.5, 5),
      coi = rbeta(ns, 2, 6)
    ),
    data = list(Y = rbeta(nobs, shape1 = 2, shape2 = 3))
  )
}

make_prep_outcome_discrete <- function(ns = 160, nobs = 12, seed = 1005) {
  set.seed(seed)
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(plogis(rnorm(ns * nobs)), ncol = nobs),
    phi = rgamma(ns, shape = 5, rate = 1),
    shape = rgamma(ns, shape = 4, rate = 1),
    sigma = rgamma(ns, shape = 3, rate = 2),
    zi = rbeta(ns, 1.5, 5)
  )
  prep$data <- list(
    Y = rpois(nobs, lambda = 5),
    trials = trials,
    lb = rep(NULL, nobs),
    ub = rep(NULL, nobs)
  )
  prep
}

expect_outcome_random <- function(family_fun, prep, i, support = c(-Inf, Inf), 
                                  check_integer = FALSE, seed = 1234) {
  set.seed(seed)
  random_outcome <- family_fun(i, prep = prep)

  testthat::expect_length(random_outcome, prep$ndraws)
  testthat::expect_true(all(is.finite(random_outcome)))

  if (is.finite(support[1])) {
    testthat::expect_true(all(random_outcome >= support[1]))
  }
  if (is.finite(support[2])) {
    testthat::expect_true(all(random_outcome <= support[2]))
  }
  if (isTRUE(check_integer)) {
    testthat::expect_true(all(abs(random_outcome - round(random_outcome)) < 1e-8))
  }
}

expect_outcome_modes <- function(family_fun, prep, i, q_ref, p_ref, 
                                 support = c(-Inf, Inf), 
                                 check_integer = FALSE) {
  expect_outcome_random(
    family_fun = family_fun, prep = prep, i = i, support = support,
    check_integer = check_integer
  )

  prob <- family_fun(i, prep = prep, output = "probability", q = q_ref)
  pit <- family_fun(i, prep = prep, output = "pit", q = q_ref)
  dens <- family_fun(i, prep = prep, output = "density", q = q_ref)
  q <- family_fun(i, prep = prep, output = "quantile", p = p_ref)

  for (x in list(prob, pit, dens, q)) {
    testthat::expect_type(x, "double")
    testthat::expect_length(x, prep$ndraws)
    testthat::expect_true(all(is.finite(x) | is.na(x)))
  }
}

test_that("posterior_predict outcome argument works for continuous families", {
  i <- 3

  family_specs <- list(
    gaussian = list(
      fun = brms:::posterior_predict_gaussian, q_ref = 0.25, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_gaussian_outcome()
    ),
    student = list(
      fun = brms:::posterior_predict_student, q_ref = 0.25, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_student_outcome()
    ),
    lognormal = list(
      fun = brms:::posterior_predict_lognormal, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    gamma = list(
      fun = brms:::posterior_predict_gamma, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    weibull = list(
      fun = brms:::posterior_predict_weibull, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    exponential = list(
      fun = brms:::posterior_predict_exponential, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    shifted_lognormal = list(
      fun = brms:::posterior_predict_shifted_lognormal, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    skew_normal = list(
      fun = brms:::posterior_predict_skew_normal, q_ref = 0.25, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_positive_outcome()
    ),
    frechet = list(
      fun = brms:::posterior_predict_frechet, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    exgaussian = list(
      fun = brms:::posterior_predict_exgaussian, q_ref = 1.2, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_positive_outcome()
    ),
    inverse_gaussian = list(
      fun = brms:::posterior_predict_inverse.gaussian, q_ref = 1.2, p_ref = 0.73,
      support = c(0, Inf), prep = make_prep_positive_outcome()
    ),
    gen_extreme_value = list(
      fun = brms:::posterior_predict_gen_extreme_value, q_ref = 0.25, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_positive_outcome()
    ),
    asym_laplace = list(
      fun = brms:::posterior_predict_asym_laplace, q_ref = 0.25, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_positive_outcome()
    ),
    zero_inflated_asym_laplace = list(
      fun = brms:::posterior_predict_zero_inflated_asym_laplace,
      q_ref = 0.25, p_ref = 0.73,
      support = c(-Inf, Inf), prep = make_prep_positive_outcome()
    ),
    xbeta = list(
      fun = brms:::posterior_predict_xbeta, q_ref = 0.4, p_ref = 0.73,
      support = c(0, 1), prep = make_prep_beta_outcome(), requires = "betareg"
    ),
    beta = list(
      fun = brms:::posterior_predict_beta, q_ref = 0.4, p_ref = 0.73,
      support = c(0, 1), prep = make_prep_beta_outcome()
    ),
    zero_inflated_beta = list(
      fun = brms:::posterior_predict_zero_inflated_beta, q_ref = 0.4, p_ref = 0.8,
      support = c(0, 1), prep = make_prep_beta_outcome()
    ),
    zero_one_inflated_beta = list(
      fun = brms:::posterior_predict_zero_one_inflated_beta, q_ref = 0.4, p_ref = 0.8,
      support = c(0, 1), prep = make_prep_beta_outcome()
    )
  )

  for (spec in family_specs) {
    if (!is.null(spec$requires)) {
      skip_if_not_installed(spec$requires)
    }
    expect_outcome_modes(
      family_fun = spec$fun,
      prep = spec$prep,
      i = i,
      q_ref = spec$q_ref,
      p_ref = spec$p_ref,
      support = spec$support
    )
  }
})

test_that("posterior_predict outcome argument works for discrete families", {
  i <- 4
  prep <- make_prep_outcome_discrete()

  family_specs <- list(
    bernoulli = list(fun = brms:::posterior_predict_bernoulli, q_ref = 1),
    binomial = list(fun = brms:::posterior_predict_binomial, q_ref = 3),
    beta_binomial = list(fun = brms:::posterior_predict_beta_binomial, q_ref = 3),
    poisson = list(fun = brms:::posterior_predict_poisson, q_ref = 3),
    negbinomial = list(fun = brms:::posterior_predict_negbinomial, q_ref = 3),
    negbinomial2 = list(fun = brms:::posterior_predict_negbinomial2, q_ref = 3),
    geometric = list(fun = brms:::posterior_predict_geometric, q_ref = 3),
    discrete_weibull = list(
      fun = brms:::posterior_predict_discrete_weibull, q_ref = 3
    ),
    com_poisson = list(fun = brms:::posterior_predict_com_poisson, q_ref = 3),
    zero_inflated_poisson = list(
      fun = brms:::posterior_predict_zero_inflated_poisson, q_ref = 3
    ),
    zero_inflated_binomial = list(
      fun = brms:::posterior_predict_zero_inflated_binomial, q_ref = 3
    ),
    zero_inflated_beta_binomial = list(
      fun = brms:::posterior_predict_zero_inflated_beta_binomial, q_ref = 3
    ),
    zero_inflated_negbinomial = list(
      fun = brms:::posterior_predict_zero_inflated_negbinomial, q_ref = 3
    )
  )

  for (spec in family_specs) {
    if (!is.null(spec$requires)) {
      skip_if_not_installed(spec$requires)
    }
    expect_outcome_modes(
      family_fun = spec$fun,
      prep = prep,
      i = i,
      q_ref = spec$q_ref,
      p_ref = 0.81,
      support = c(0, Inf),
      check_integer = TRUE
    )
  }
})

test_that("compute_cdf returns correct CDF for non-truncated distributions", {
  # Non-truncated, non-randomized: raw CDF F(q)
  q <- 3
  out <- brms:::compute_cdf(q = q, dist = "pois", lb = NULL, ub = NULL,
  randomized = FALSE, lambda = 5)
  expect_equal(out, ppois(q, lambda = 5))

  q <- 2
  out <- brms:::compute_cdf(q = q, dist = "binom", lb = NULL, ub = NULL,
  randomized = FALSE, size = 10, prob = 0.5)
  expect_equal(out, pbinom(q, size = 10, prob = 0.5))
})

test_that("compute_cdf with randomized = TRUE returns value in [F(q-1), F(q)]", {
  # Randomized PIT: F(y-1) + V * [F(y) - F(y-1)] with V ~ Unif(0,1)
  set.seed(42)
  q <- 5
  Fq <- ppois(q, lambda = 3)
  Fqm1 <- ppois(q - 1, lambda = 3)

  out <- brms:::compute_cdf(q = q, dist = "pois", lb = NULL, ub = NULL, 
                            randomized = TRUE, lambda = 3)
  expect_true(out >= Fqm1)
  expect_true(out <= Fq)
})

test_that("compute_cdf with randomized = TRUE and truncation returns value in 
valid range", {
  set.seed(123)
  q <- 4
  lb <- 2
  ub <- 7
  denom <- ppois(ub, lambda = 5) - ppois(lb, lambda = 5)
  Fq <- (ppois(q, lambda = 5) - ppois(lb, lambda = 5)) / denom
  Fqm1 <- (ppois(q - 1, lambda = 5) - ppois(lb, lambda = 5)) / denom

  out <- brms:::compute_cdf(q = q, dist = "pois", lb = lb, ub = ub, 
                            randomized = TRUE, lambda = 5)
  expect_true(out >= Fqm1)
  expect_true(out <= Fq)
  expect_true(out >= 0 && out <= 1)
})

test_that("compute_cdf handles zero denominator (lb == ub) without unexpected 
behaviour", {
  q <- 3
  lb <- 1
  ub <- 1

  out <- tryCatch(
    brms:::compute_cdf(q = q, dist = "pois", lb = lb, ub = ub, 
                       randomized = FALSE, lambda = 2),
    error = function(e) structure(list(error = TRUE, message = e$message))
  )

  if (is.list(out) && isTRUE(out$error)) {
    expect_true(grepl("zero|denom|divide|trunc", out$message, ignore.case = TRUE))
  } else {
    expect_true(is.nan(out) || is.na(out))
  }
})


test_that("compute_cdf respects lower.tail and log.p", {
  q <- 0.4
  mu <- 0.1
  sd <- 1.7

  base <- brms:::compute_cdf(
    q = q, distribution = "norm", lb = NULL, ub = NULL, randomized = FALSE,
    lower.tail = TRUE, log.p = FALSE,
    mean = mu, sd = sd
  )

  upper <- brms:::compute_cdf(
    q = q, distribution = "norm", lb = NULL, ub = NULL, randomized = FALSE,
    lower.tail = FALSE, log.p = FALSE,
    mean = mu, sd = sd
  )

  log_base <- brms:::compute_cdf(
    q = q, distribution = "norm", lb = NULL, ub = NULL, randomized = FALSE,
    lower.tail = TRUE, log.p = TRUE,
    mean = mu, sd = sd
  )

  log_upper <- brms:::compute_cdf(
    q = q, distribution = "norm", lb = NULL, ub = NULL, randomized = FALSE,
    lower.tail = FALSE, log.p = TRUE,
    mean = mu, sd = sd
  )

  expected <- pnorm(q, mean = mu, sd = sd)

  expect_equal(base, expected)
  expect_equal(upper, 1 - expected)
  expect_equal(log_base, log(expected))
  expect_equal(log_upper, log(1 - expected))
})

test_that("posterior_predict forwards lower.tail and log.p correctly", {
  fit <- brms:::rename_pars(brms:::brmsfit_example3)
  q_ref <- 15

  p_lower <- posterior_predict(fit, output = "probability", q = q_ref,
    lower.tail = TRUE, log.p = FALSE)

  p_upper <- posterior_predict(fit, output = "probability", q = q_ref,
    lower.tail = FALSE, log.p = FALSE)

  log_lower <- posterior_predict(fit, output = "probability", q = q_ref,
    lower.tail = TRUE, log.p = TRUE)

  log_upper <- posterior_predict(fit, output = "probability", q = q_ref,
    lower.tail = FALSE, log.p = TRUE)
  
  expect_equal(dim(p_lower), dim(p_upper))
  expect_equal(p_upper, 1 - p_lower)
  expect_equal(log_lower, log(p_lower))
  expect_equal(log_upper, log(p_upper))
})