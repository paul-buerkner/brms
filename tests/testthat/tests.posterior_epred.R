context("Tests for posterior_epred helper functions")

# to reduce testing time on CRAN
skip_on_cran()

test_that("Cox RMST controls flow through preparation and fitted wrappers", {
  skip_if_not_installed("splines2")
  skip_if_not_installed("statmod")
  fit <- brm(y ~ 1, data = data.frame(y = 1:10), family = cox(), empty = TRUE)
  fit$fit <- brms:::rename_pars(brms:::brmsfit_example1)$fit
  for (i in 1:5) {
    fit <- brms:::add_dummy_draws(
      fit, paste0("sbhaz[", i, "]"), mean = 0.2, sd = 0
    )
  }
  prep <- prepare_predictions(
    fit, draw_ids = 1:3, horizon = 2, nodes = 16, check_response = FALSE
  )
  expect_equal(prep$bhaz$horizon, 2)
  expect_equal(prep$bhaz$nodes, 16)
  expected <- posterior_epred(prep, dpar = NULL, nlpar = NULL, sort = FALSE)
  expect_equal(
    posterior_epred(fit, draw_ids = 1:3, horizon = 2, nodes = 16), expected
  )
  expect_equal(
    fitted(fit, draw_ids = 1:3, horizon = 2, nodes = 16),
    posterior_summary(expected)
  )
  expect_s3_class(pp_check(
    fit, type = "error_binned", x = "y", draw_ids = 1:3,
    horizon = 2, nodes = 16
  ), "ggplot")

  prep <- prepare_predictions(fit, draw_ids = 1:3, check_response = FALSE)
  expect_true(all(is.finite(posterior_predict(prep))))
  expect_true(all(is.finite(posterior_epred(
    fit, draw_ids = 1:3, dpar = "mu"
  ))))
  expect_true(all(is.finite(fitted(fit, draw_ids = 1:3, scale = "linear"))))
  expect_error(
    posterior_epred(prep, dpar = NULL, nlpar = NULL, sort = FALSE), "horizon"
  )
  other <- brms:::rename_pars(brms:::brmsfit_example1)
  expect_equal(
    posterior_epred(other, draw_ids = 1:3, horizon = 2, nodes = 16),
    posterior_epred(other, draw_ids = 1:3)
  )
})

test_that("Prepared expected predictions reject late preparation arguments", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep()
  for (args in list(list(horizon = 2), list(nodes = 16), list(newdata = NULL))) {
    expect_error(
      do.call(posterior_epred, c(
        list(prep, dpar = NULL, nlpar = NULL, sort = FALSE), args
      )),
      "must be supplied when preparing"
    )
  }
})

test_that("Cox RMST agrees with independent integration and truncation", {
  skip_if_not_installed("splines2")
  skip_if_not_installed("statmod")
  prep <- cox_test_prep(ns = 3L, nobs = 1L)
  sb <- prep$bhaz$sbhaz[[1]][1, , drop = FALSE]
  basis <- prep$bhaz$basis
  for (horizon in c(0, 2, 6)) {
    prep$bhaz$horizon <- horizon
    got <- posterior_epred(prep, dpar = NULL, nlpar = NULL, sort = FALSE)
    for (i in seq_len(3)) {
      mu <- prep$dpars$mu[i, 1]
      reference <- integrate(function(t) {
        H <- as.vector(brms:::bhaz_basis_matrix(
          pmin(t, 4), basis = basis, integrate = TRUE
        ) %*% as.vector(sb))
        rate <- as.vector(brms:::bhaz_basis_matrix(4, basis = basis) %*%
                            as.vector(sb))
        exp(-mu * (H + pmax(t - 4, 0) * rate))
      }, 0, horizon, rel.tol = 1e-10, subdivisions = 1000)$value
      expect_equal(got[i, 1], reference, tolerance = 1e-7)
    }
  }
  prep$data$lb <- 1
  prep$data$ub <- 2
  prep$bhaz$horizon <- 0.5
  expect_equal(as.vector(posterior_epred(
    prep, dpar = NULL, nlpar = NULL, sort = FALSE
  )), rep(0.5, 3))
  set.seed(42)
  pred <- replicate(2000, posterior_predict(prep)[1, 1])
  prep$bhaz$horizon <- 1.5
  got <- posterior_epred(prep, dpar = NULL, nlpar = NULL, sort = FALSE)
  expect_equal(got[1, 1], mean(pmin(pred, 1.5)), tolerance = 0.02)
  prep$bhaz$horizon <- NULL
  expect_error(posterior_epred(prep, dpar = NULL, nlpar = NULL,
                               sort = FALSE), "horizon")
  prep$bhaz$tail <- "none"
  prep$bhaz$horizon <- 5
  expect_error(posterior_epred(prep, dpar = NULL,
                               nlpar = NULL, sort = FALSE), "support")
  prep <- cox_test_prep(ns = 2L)
  prep$data$lb <- c(0.25, 1, 2)
  prep$data$ub <- c(2.5, 3.5, Inf)
  prep$bhaz$tail <- "loglinear"
  prep$bhaz$horizon <- 6
  got <- posterior_epred(prep, dpar = NULL, nlpar = NULL, sort = FALSE)
  for (i in seq_len(prep$nobs)) {
    ref <- integrate(function(t) brms:::pcox(
      t, prep$dpars$mu[1, i], sbhaz = prep$bhaz$sbhaz[[1]][1, , drop = FALSE],
      basis = prep$bhaz$basis, tail = "loglinear", lower.tail = FALSE,
      lb = prep$data$lb[i], ub = prep$data$ub[i]
    ), 0, 6, subdivisions = 1000, rel.tol = 1e-9)$value
    expect_equal(got[1, i], ref, tolerance = 1e-7)
  }
  prep$dpars$mu[] <- -1
  expect_error(posterior_epred(prep, dpar = NULL, nlpar = NULL,
                               sort = FALSE), "finite positive")
})

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

test_that("truncated discrete posterior_epred keeps the lower bound", {
  # the mean is taken over lb:ub and normalized by F(ub) - F(lb - 1), as in
  # the generated Stan code and in log_lik; see #1923
  epred_prep <- function(mu, family, ...) {
    prep <- structure(list(ndraws = 2L, nobs = 1L), class = "brmsprep")
    prep$dpars <- c(list(mu = matrix(mu, nrow = 2, ncol = 1)), list(...))
    prep$family <- family
    prep
  }
  bnd <- function(x) matrix(x, nrow = 2, ncol = 1)
  lambda <- 3

  got <- suppressMessages(brms:::posterior_epred_trunc_poisson(
    epred_prep(lambda, poisson()), lb = bnd(2), ub = bnd(6)
  ))
  expect_equal(
    got[1],
    sum((2:6) * dpois(2:6, lambda)) /
      (ppois(6, lambda) - ppois(1, lambda))
  )

  # a lower bound at the edge of the support sends lb - 1 outside it
  got <- suppressMessages(brms:::posterior_epred_trunc_poisson(
    epred_prep(lambda, poisson()), lb = bnd(0), ub = bnd(6)
  ))
  expect_equal(
    got[1],
    sum((0:6) * dpois(0:6, lambda)) / ppois(6, lambda)
  )

  # a non-integer bound rounds to the support it admits, as in log_lik.
  # min_lb is a grid origin shared by all observations, so one fractional
  # bound would otherwise zero every column
  got <- suppressMessages(brms:::posterior_epred_trunc_poisson(
    epred_prep(lambda, poisson()), lb = bnd(2.5), ub = bnd(6)
  ))
  expect_equal(
    got[1],
    sum((3:6) * dpois(3:6, lambda)) /
      (ppois(6, lambda) - ppois(2, lambda))
  )

  # negbinomial covers the other cdf shape in the affected set
  got <- suppressMessages(brms:::posterior_epred_trunc_negbinomial(
    epred_prep(lambda, negbinomial(), shape = matrix(2, nrow = 2, ncol = 1)),
    lb = bnd(2), ub = bnd(6)
  ))
  expect_equal(
    got[1],
    sum((2:6) * dnbinom(2:6, size = 2, mu = lambda)) /
      (pnbinom(6, size = 2, mu = lambda) - pnbinom(1, size = 2, mu = lambda))
  )
})
