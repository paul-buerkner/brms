context("Validation tests for new discrete quantile distributions")

test_that("qbeta_binomial satisfies the discrete quantile definition", {
  skip_if_not_installed("extraDistr")

  size <- 12
  mu <- 0.35
  phi <- 8
  p <- c(0.01, 0.1, 0.5, 0.85, 0.99)

  q <- brms:::qbeta_binomial(p, size = size, mu = mu, phi = phi)
  F_q <- brms:::pbeta_binomial(q, size = size, mu = mu, phi = phi)
  F_qm1 <- brms:::pbeta_binomial(q - 1, size = size, mu = mu, phi = phi)

  expect_true(all(F_q >= p))
  expect_true(all(F_qm1 < p))
  expect_true(all(q >= 0 & q <= size))
})

test_that("qcom_poisson satisfies the discrete quantile definition", {
  mu <- 2
  shape <- 0.8
  p <- c(0.01, 0.1, 0.5, 0.85, 0.99)

  q <- brms:::qcom_poisson(p, mu = mu, shape = shape)
  F_q <- brms:::pcom_poisson(q, mu = mu, shape = shape)
  F_qm1 <- brms:::pcom_poisson(q - 1, mu = mu, shape = shape)

  expect_true(all(F_q >= p))
  expect_true(all(ifelse(q > 0, F_qm1 < p, TRUE)))
  expect_true(all(q >= 0))
})

test_that("qcom_poisson reduces to qpois when shape is 1", {
  mu <- 3
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    brms:::qcom_poisson(p, mu = mu, shape = 1),
    qpois(p, lambda = mu)
  )
})

test_that("qhurdle_poisson satisfies the discrete quantile definition", {
  lambda <- 2
  hu <- 0.2
  p <- c(0.01, 0.1, 0.5, 0.85, 0.99)

  q <- brms:::qhurdle_poisson(p, lambda = lambda, hu = hu)
  F_q <- brms:::phurdle_poisson(q, lambda = lambda, hu = hu)
  F_qm1 <- brms:::phurdle_poisson(q - 1, lambda = lambda, hu = hu)

  expect_true(all(F_q >= p))
  expect_true(all(ifelse(q > 0, F_qm1 < p, TRUE)))
  expect_true(all(q >= 0))
})

test_that("qzero_inflated_negbinomial satisfies the discrete quantile definition", {
  mu <- 4
  shape <- 2.5
  zi <- 0.3
  p <- c(0.01, 0.1, 0.5, 0.85, 0.99)

  q <- brms:::qzero_inflated_negbinomial(p, mu = mu, shape = shape, zi = zi)
  F_q <- brms:::pzero_inflated_negbinomial(q, mu = mu, shape = shape, zi = zi)
  F_qm1 <- brms:::pzero_inflated_negbinomial(q - 1, mu = mu, shape = shape, zi = zi)

  expect_true(all(F_q >= p))
  expect_true(all(F_qm1 < p))
  expect_true(all(q >= 0))
})

test_that("quantile functions are monotone in probability", {
  p <- seq(0.01, 0.99, length.out = 50)

  q_bb <- brms:::qbeta_binomial(p, size = 20, mu = 0.45, phi = 6)
  q_zinb <- brms:::qzero_inflated_negbinomial(p, mu = 5, shape = 2, zi = 0.25)
  q_comp <- brms:::qcom_poisson(p, mu = 2, shape = 1.5)

  expect_true(all(diff(q_bb) >= 0))
  expect_true(all(diff(q_zinb) >= 0))
  expect_true(all(diff(q_comp) >= 0))
})

test_that("quantile functions respect lower.tail and log.p", {
  p_ref <- 0.2

  q_bb_upper <- brms:::qbeta_binomial(
    p_ref, size = 20, mu = 0.45, phi = 6,
    lower.tail = FALSE, log.p = FALSE
  )
  q_bb_ref <- brms:::qbeta_binomial(
    1 - p_ref, size = 20, mu = 0.45, phi = 6,
    lower.tail = TRUE, log.p = FALSE
  )
  expect_equal(q_bb_upper, q_bb_ref)

  q_zinb_log <- brms:::qzero_inflated_negbinomial(
    log(p_ref), mu = 5, shape = 2, zi = 0.25,
    lower.tail = TRUE, log.p = TRUE
  )
  q_zinb_ref <- brms:::qzero_inflated_negbinomial(
    p_ref, mu = 5, shape = 2, zi = 0.25,
    lower.tail = TRUE, log.p = FALSE
  )
  expect_equal(q_zinb_log, q_zinb_ref)

  q_comp_upper <- brms:::qcom_poisson(
    p_ref, mu = 2, shape = 0.8,
    lower.tail = FALSE, log.p = FALSE
  )
  q_comp_lower <- brms:::qcom_poisson(
    1 - p_ref, mu = 2, shape = 0.8,
    lower.tail = TRUE, log.p = FALSE
  )
  expect_equal(q_comp_upper, q_comp_lower)

  q_comp_log <- brms:::qcom_poisson(
    log(p_ref), mu = 2, shape = 0.8,
    lower.tail = TRUE, log.p = TRUE
  )
  q_comp_ref <- brms:::qcom_poisson(
    p_ref, mu = 2, shape = 0.8,
    lower.tail = TRUE, log.p = FALSE
  )
  expect_equal(q_comp_log, q_comp_ref)
})
