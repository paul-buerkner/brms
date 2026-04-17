context("Validation tests for posterior_predict density and quantile outputs")

skip_if_density_quantile_not_available <- function() {
  ns <- asNamespace("brms")
  has_density <- exists("compute_density", envir = ns, inherits = FALSE)
  has_quantile <- exists("compute_quantile", envir = ns, inherits = FALSE)
  if (!(has_density && has_quantile)) {
    testthat::skip("Density/quantile posterior_predict helpers are not available.")
  }
}

make_prep_gaussian <- function(ns = 120, nobs = 8) {
  set.seed(1001)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rgamma(ns, shape = 4, rate = 3)
  )
  prep$data <- list(
    Y = rnorm(nobs),
    lb = rep(NULL, nobs),
    ub = rep(NULL, nobs)
  )
  prep
}

make_prep_student <- function(ns = 120, nobs = 8) {
  set.seed(1002)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rgamma(ns, shape = 4, rate = 3),
    nu = rgamma(ns, shape = 6, rate = 1) + 2
  )
  prep$data <- list(
    Y = rnorm(nobs),
    lb = rep(NULL, nobs),
    ub = rep(NULL, nobs)
  )
  prep
}

make_prep_count <- function(ns = 150, nobs = 10) {
  set.seed(1003)
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

test_that("continuous families: density matches finite-difference CDF", {

  h <- 1e-4
  i <- 3
  q_ref <- 0.25

  prep_g <- make_prep_gaussian()
  f_g <- brms:::posterior_predict_gaussian(i, 
    prep = prep_g, output = "density", q = q_ref)
  F_plus_g <- brms:::posterior_predict_gaussian(i, 
    prep = prep_g, output = "probability", q = q_ref + h)
  F_minus_g <- brms:::posterior_predict_gaussian(i, 
    prep = prep_g, output = "probability", q = q_ref - h)
  fd_g <- (F_plus_g - F_minus_g) / (2 * h)
  expect_equal(f_g, fd_g, tolerance = 1e-3)

  prep_t <- make_prep_student()
  f_t <- brms:::posterior_predict_student(i, 
    prep = prep_t, output = "density", q = q_ref)
  F_plus_t <- brms:::posterior_predict_student(i, 
    prep = prep_t, output = "probability", q = q_ref + h)
  F_minus_t <- brms:::posterior_predict_student(i, 
    prep = prep_t, output = "probability", q = q_ref - h)
  fd_t <- (F_plus_t - F_minus_t) / (2 * h)
  expect_equal(f_t, fd_t, tolerance = 2e-3)
})

test_that("density output supports log = TRUE", {

  i <- 3
  q_ref <- 0.25
  prep_g <- make_prep_gaussian()

  dens <- brms:::posterior_predict_gaussian(
    i, prep = prep_g, output = "density", q = q_ref, log = FALSE
  )
  log_dens <- brms:::posterior_predict_gaussian(
    i, prep = prep_g, output = "density", q = q_ref, log = TRUE
  )

  expect_equal(log_dens, log(dens), tolerance = 1e-10)
})

test_that("discrete families: density matches CDF increments", {

  i <- 4
  q_ref <- 3
  prep <- make_prep_count()

  check_density_increment <- function(family_fun, ...) {
    dens <- family_fun(i, prep = prep, output = "density", q = q_ref, ...)
    F_q <- family_fun(i, prep = prep, output = "probability", q = q_ref, ...)
    F_qm1 <- family_fun(i, prep = prep, output = "probability", q = q_ref - 1, ...)
    expect_equal(dens, F_q - F_qm1, tolerance = 1e-10)
  }

  check_density_increment(brms:::posterior_predict_binomial)
  check_density_increment(brms:::posterior_predict_beta_binomial)
  check_density_increment(brms:::posterior_predict_poisson)
  check_density_increment(brms:::posterior_predict_negbinomial)
  check_density_increment(brms:::posterior_predict_negbinomial2)
  check_density_increment(brms:::posterior_predict_zero_inflated_negbinomial)
})

test_that("quantile output inverts probability output (continuous)", {

  i <- 2
  p_ref <- 0.73

  prep_g <- make_prep_gaussian()
  q_g <- brms:::posterior_predict_gaussian(i, prep = prep_g, 
    output = "quantile", p = p_ref)
  F_q_g <- brms:::posterior_predict_gaussian(i, prep = prep_g, 
    output = "probability", q = q_g)
  expect_equal(F_q_g, rep(p_ref, length(F_q_g)), tolerance = 1e-8)

  prep_t <- make_prep_student()
  q_t <- brms:::posterior_predict_student(i, prep = prep_t, 
    output = "quantile", p = p_ref)
  F_q_t <- brms:::posterior_predict_student(i, prep = prep_t, 
    output = "probability", q = q_t)
  expect_equal(F_q_t, rep(p_ref, length(F_q_t)), tolerance = 1e-8)
})

test_that("quantile output inverts CDF inequalities (discrete)", {
  skip_if_density_quantile_not_available()

  i <- 5
  p_ref <- 0.81
  prep <- make_prep_count()

  check_quantile_discrete <- function(family_fun, ...) {
    q <- family_fun(i, prep = prep, output = "quantile", p = p_ref, ...)
    F_q <- family_fun(i, prep = prep, output = "probability", q = q, ...)
    F_qm1 <- family_fun(i, prep = prep, output = "probability", q = q - 1, ...)

    expect_true(all(F_q >= p_ref))
    expect_true(all(F_qm1 < p_ref))
  }

  check_quantile_discrete(brms:::posterior_predict_binomial)
  check_quantile_discrete(brms:::posterior_predict_beta_binomial)
  check_quantile_discrete(brms:::posterior_predict_poisson)
  check_quantile_discrete(brms:::posterior_predict_negbinomial)
  check_quantile_discrete(brms:::posterior_predict_negbinomial2)
  check_quantile_discrete(brms:::posterior_predict_zero_inflated_negbinomial)
})

test_that("quantile output supports lower.tail and log.p", {

  i <- 2
  prep_g <- make_prep_gaussian()
  p_ref <- 0.2

  q_upper <- brms:::posterior_predict_gaussian(
    i, prep = prep_g, output = "quantile",
    p = p_ref, lower.tail = FALSE, log.p = FALSE
  )
  F_upper <- brms:::posterior_predict_gaussian(
    i, prep = prep_g, output = "probability",
    q = q_upper, lower.tail = FALSE, log.p = FALSE
  )
  expect_equal(F_upper, rep(p_ref, length(F_upper)), tolerance = 1e-8)

  q_log <- brms:::posterior_predict_gaussian(
    i, prep = prep_g, output = "quantile",
    p = log(p_ref), lower.tail = TRUE, log.p = TRUE
  )
  F_log <- brms:::posterior_predict_gaussian(
    i, prep = prep_g, output = "probability",
    q = q_log, lower.tail = TRUE, log.p = FALSE
  )
  expect_equal(F_log, rep(p_ref, length(F_log)), tolerance = 1e-8)
})

test_that("random baseline agrees with density for selected discrete families", {

  set.seed(1234)
  i <- 6
  q_ref <- 2
  prep <- make_prep_count(ns = 4000, nobs = 10)

  check_random_baseline <- function(family_fun, ...) {
    draws <- family_fun(i, prep = prep, output = "random", ...)
    p_emp <- mean(draws == q_ref)
    p_dens <- mean(family_fun(i, prep = prep, output = "density", q = q_ref, ...))
    expect_equal(p_emp, p_dens, tolerance = 0.03)
  }

  check_random_baseline(brms:::posterior_predict_binomial)
  check_random_baseline(brms:::posterior_predict_poisson)
  check_random_baseline(brms:::posterior_predict_negbinomial)
})
