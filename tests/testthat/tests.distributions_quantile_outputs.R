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

  expect_true(all(diff(q_bb) >= 0))
  expect_true(all(diff(q_zinb) >= 0))
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
})
