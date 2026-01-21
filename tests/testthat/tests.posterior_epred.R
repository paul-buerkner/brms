context("Tests for posterior_epred helper functions")

# to reduce testing time on CRAN
skip_on_cran()

test_that("posterior_epred helper functions run without errors", {
  # actually run posterior_epred.brmsfit that call the helper functions
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  add_dummy_draws <- brms:::add_dummy_draws
  fit <- add_dummy_draws(fit, "shape", dist = "exp")
  fit <- add_dummy_draws(fit, "alpha", dist = "norm")
  fit <- add_dummy_draws(fit, "hu", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_dummy_draws(fit, "phi", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_dummy_draws(fit, "zi", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_dummy_draws(fit, "quantile", dist = "beta", shape1 = 2, shape2 = 1)
  fit <- add_dummy_draws(fit, "xi", dist = "unif", min = -1, max = 0.5)
  fit <- add_dummy_draws(fit, "ndt", dist = "exp")
  fit$formula$formula <- update(fit$formula$formula, .~. - arma(visit, patient))
  prep <- brms:::prepare_predictions(fit)
  prep$dpars$mu <- brms:::get_dpar(prep, "mu")
  prep$dpars$sigma <- brms:::get_dpar(prep, "sigma")
  prep$dpars$nu <- brms:::get_dpar(prep, "nu")
  ndraws <- ndraws(fit)
  nobs <- nobs(fit)

  # test preparation of truncated models
  prep$data$lb <- 0
  prep$data$ub <- 200
  mu <- brms:::posterior_epred_trunc(prep)
  expect_equal(dim(mu), c(ndraws, nobs))

  # pseudo log-normal model
  fit$family <- fit$formula$family <- lognormal()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)),
               c(ndraws, nobs))

  # pseudo shifted log-normal model
  fit$family <- fit$formula$family <- shifted_lognormal()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)),
               c(ndraws, nobs))

  # pseudo skew-normal model
  fit$family <- fit$formula$family <- skew_normal()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)),
               c(ndraws, nobs))

  # pseudo asym_laplace model
  fit$family <- fit$formula$family <- asym_laplace()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)),
               c(ndraws, nobs))

  # pseudo zero_inflated_asym_laplace model
  fit$family <- fit$formula$family <- brmsfamily("zero_inflated_asym_laplace")
  expect_equal(dim(posterior_epred(fit, summary = FALSE)),
               c(ndraws, nobs))

  # pseudo gen_extreme_value model
  fit$family <- fit$formula$family <- gen_extreme_value()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)),
               c(ndraws, nobs))

  # pseudo weibull model
  fit$formula$pforms <- NULL
  fit$family <- fit$formula$family <- weibull()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(ndraws, nobs))

  # pseudo binomial model
  old_formula <- fit$formula$formula
  fit$formula$formula <- update(fit$formula$formula, . | trials(100) ~ .)
  fit$autocor <- brms:::cor_empty()
  fit$family <- fit$formula$family <- binomial()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(ndraws, nobs))

  # pseudo beta-binomial model
  fit$family <- fit$formula$family <- beta_binomial()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(ndraws, nobs))

  # pseudo zero inflated binomial model
  fit$family <- fit$formula$family <- zero_inflated_binomial()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(ndraws, nobs))

  # pseudo zero inflated beta binomial model
  fit$family <- fit$formula$family <- zero_inflated_beta_binomial()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(ndraws, nobs))

  # pseudo hurdle poisson model
  fit$formula$formula <- old_formula
  fit$family <- fit$formula$family <- hurdle_poisson()
  fit$formula <- bf(count ~ Trt*Age + mo(Exp) + offset(Age) + (1+Trt|visit),
                    family = family(fit))
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), c(ndraws, nobs))

  # pseudo zero-inflated poisson model
  fit$family <- fit$formula$family <- zero_inflated_poisson()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), c(ndraws, nobs))

  # pseudo custom model
  posterior_epred_test <- function(prep) {
    prep$dpars$mu
  }
  fit$family <- fit$formula$family <- custom_family(
    "test", dpars = "mu", links = c("logit"),
    type = "int", vars = "trials[n]"
  )
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), c(ndraws, nobs))

  # truncated continuous models
  prep$dpars$shape <- c(as.matrix(fit, variable = "shape"))
  mu <- brms:::posterior_epred_trunc_gaussian(prep, lb = 0, ub = 10)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- brms:::posterior_epred_trunc_student(prep, lb = -Inf, ub = 15)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- brms:::posterior_epred_trunc_lognormal(prep, lb = 2, ub = 15)
  expect_equal(dim(mu), c(ndraws, nobs))

  prep$dpars$mu <- exp(prep$dpars$mu)
  mu <- brms:::posterior_epred_trunc_gamma(prep, lb = 1, ub = 7)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- brms:::posterior_epred_trunc_exponential(prep, lb = 0, ub = Inf)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- SW(brms:::posterior_epred_trunc_weibull(prep, lb = -Inf, ub = Inf))
  expect_equal(dim(mu), c(ndraws, nobs))

  # truncated discrete models
  data <- list(Y = sample(100, 10), trials = 1:10, N = 10)
  lb <- matrix(0, nrow = ndraws, ncol = nobs)
  ub <- matrix(100, nrow = ndraws, ncol = nobs)
  mu <- brms:::posterior_epred_trunc_poisson(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- brms:::posterior_epred_trunc_negbinomial(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- brms:::posterior_epred_trunc_negbinomial2(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(ndraws, nobs))

  mu <- brms:::posterior_epred_trunc_geometric(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(ndraws, nobs))

  prep$data$trials <- 120
  lb <- matrix(-Inf, nrow = ndraws, ncol = nobs)
  prep$dpars$mu <- brms:::inv_link(prep$dpars$mu, "logit")
  mu <- brms:::posterior_epred_trunc_binomial(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(ndraws, nobs))
})

test_that("posterior_epred_lagsar runs without errors", {
  prep <- list(
    dpars = list(mu = matrix(rnorm(30), nrow = 3)),
    ac = list(
      lagsar = matrix(c(0.3, 0.5, 0.7)),
      Msar = matrix(1:100, 10, 10)
    ),
    ndraws = 3,
    nobs = 10,
    family = gaussian()
  )
  mu_new <- brms:::posterior_epred_lagsar(prep)
  expect_equal(dim(mu_new), dim(prep$dpars$mu))
  expect_true(!identical(mu_new, prep$dpars$mu))
})

test_that("posterior_epred for advanced count data distributions runs without errors", {
  ns <- 15
  nobs <- 5
  ncat <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = array(rbeta(ns*nobs, 2, 2), dim = c(ns, nobs)),
    shape = array(rexp(ns*nobs, 3), dim = c(ns, nobs))
  )
  prep$family <- brmsfamily("discrete_weibull")
  pred <- suppressWarnings(brms:::posterior_epred_discrete_weibull(prep))
  expect_equal(dim(pred), c(ns, nobs))

  prep$family <- brmsfamily("com_poisson")
  pred <- suppressWarnings(brms:::posterior_epred_com_poisson(prep))
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_epred for multinomial, dirichlet_multinomial and dirichlet models runs without errors", {
  ns <- 15
  nobs <- 8
  ncat <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu1 = array(rnorm(ns*nobs), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs), dim = c(ns, nobs))
  )
  prep$data <- list(ncat = ncat, trials = sample(1:20, nobs))
  prep$refcat <- 1

  prep$family <- multinomial()
  pred <- brms:::posterior_epred_multinomial(prep = prep)
  expect_equal(dim(pred), c(ns, nobs, ncat))

  prep$family <- dirichlet_multinomial()
  pred <- brms:::posterior_epred_dirichlet_multinomial(prep = prep)
  expect_equal(dim(pred), c(ns, nobs, ncat))

  prep$family <- dirichlet()
  pred <- brms:::posterior_epred_dirichlet(prep = prep)
  expect_equal(dim(pred), c(ns, nobs, ncat))

  prep$family <- brmsfamily("dirichlet2")
  prep$dpars$mu1 <- array(rexp(ns*nobs, 1), dim = c(ns, nobs))
  prep$dpars$mu2 <- array(rexp(ns*nobs, 1), dim = c(ns, nobs))
  prep$dpars$mu3 <- array(rexp(ns*nobs, 1), dim = c(ns, nobs))
  pred <- brms:::posterior_epred_dirichlet2(prep = prep)
  expect_equal(dim(pred), c(ns, nobs, ncat))
})

test_that("posterior_epred_xbeta runs without errors", {
  ns <- 50
  nobs <- 8
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rbeta(ns * nobs, 1.2, 2.3), ncol = nobs),
    phi = rexp(ns, 0.01),
    kappa = rexp(ns, 2)
  )
  prep$data <- list(Y = rbeta(nobs, 2, 3))
  mu_new <- brms:::posterior_epred_xbeta(prep)
  expect_equal(dim(mu_new), dim(prep$dpars$mu))
  expect_true(!identical(mu_new, prep$dpars$mu))
})

test_that("posterior_epred_ordbeta runs without errors", {
  ns <- 50
  nobs <- 8
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    phi = matrix(rexp(ns * nobs, 0.1), ncol = nobs),
    cutzero = matrix(rnorm(ns * nobs, -1, 0.5), ncol = nobs),
    cutone = matrix(rnorm(ns * nobs, 0, 0.5), ncol = nobs)
  )
  prep$data <- list(Y = c(0, 0.3, 0.5, 0.7, 1, 0.2, 0.8, 0.4))
  mu_new <- brms:::posterior_epred_ordbeta(prep)
  expect_equal(dim(mu_new), dim(prep$dpars$mu))
  expect_true(all(mu_new >= 0 & mu_new <= 1))
})

test_that("posterior_epred can be reproduced by using d<family>()", {
  fit4 <- rename_pars(brms:::brmsfit_example4)
  epred4 <- posterior_epred(fit4)

  eta4 <- posterior_linpred(fit4)
  bprep4 <- prepare_predictions(fit4)
  thres4 <- bprep4$thres$thres
  disc4 <- bprep4$dpars$disc$fe$b %*% t(bprep4$dpars$disc$fe$X)
  disc4 <- exp(disc4)
  epred4_ch <- aperm(sapply(seq_len(dim(eta4)[2]), function(i) {
    dsratio(seq_len(ncol(thres4) + 1), eta4[, i, ], thres4, disc4[, i])
  }, simplify = "array"), perm = c(1, 3, 2))

  expect_equivalent(epred4, epred4_ch)
})

