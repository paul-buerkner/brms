context("Tests for numeric helper functions")

test_that("log1m_exp is accurate for arguments close to zero", {
  # the naive log1p(-exp(x)) returns -Inf for all of these, because exp(x)
  # rounds to exactly 1 and log1p(-1) is -Inf
  x <- c(-1e-18, -1e-40, -5.489115e-85)
  expect_true(all(is.infinite(log1p(-exp(x)))))
  # for |x| this small, log(1 - exp(x)) equals log(-x) to full double precision
  expect_equal(brms:::log1m_exp(x), log(-x))
})

test_that("log1m_exp matches the naive form where the naive form is valid", {
  x <- c(-0.05, -0.5, -1, -5, -30)
  expect_equal(brms:::log1m_exp(x), log(1 - exp(x)))
})

test_that("log1m_exp follows the domain conventions of Stan", {
  # Stan's log1m_exp returns NaN for arguments >= 0
  expect_true(is.nan(brms:::log1m_exp(0)))
  expect_true(is.nan(brms:::log1m_exp(1)))
  expect_true(is.na(brms:::log1m_exp(NA)))
  # the shape of the input is preserved
  m <- matrix(c(-0.1, -2, -0.5, -3), nrow = 2)
  expect_equal(dim(brms:::log1m_exp(m)), c(2L, 2L))
})

test_that("log_diff_exp survives underflow of both arguments", {
  # exp() of either argument underflows to 0, so log(exp(x) - exp(y)) is -Inf
  expect_equal(log(exp(-800) - exp(-900)), -Inf)
  expect_equal(brms:::log_diff_exp(-800, -900), -800)
  expect_equal(brms:::log_diff_exp(-1e5, -2e5), -1e5)
  # a difference that is resolvable is still resolved
  expect_equal(brms:::log_diff_exp(-800, -800.5), -800 + log1p(-exp(-0.5)))
})

test_that("log_diff_exp matches the naive form where the naive form is valid", {
  x <- c(2, 0, -3, 10)
  y <- c(1, -1, -4, 9)
  expect_equal(brms:::log_diff_exp(x, y), log(exp(x) - exp(y)))
})

test_that("log_diff_exp follows the edge cases of Stan", {
  # equal arguments give log(0) rather than NaN
  expect_equal(brms:::log_diff_exp(3, 3), -Inf)
  expect_equal(brms:::log_diff_exp(-Inf, -Inf), -Inf)
  # a negative difference is undefined
  expect_true(is.nan(brms:::log_diff_exp(1, 2)))
  expect_true(is.nan(brms:::log_diff_exp(Inf, Inf)))
  expect_equal(brms:::log_diff_exp(0, -Inf), 0)
})

test_that("log_inv_logit and log1m_inv_logit stay exact in both tails", {
  # inv_logit(x) rounds to 1 from x = 36.74 and to 0 below x = -709.78, so
  # taking the log afterwards gave exactly 0 or -Inf. The reference never
  # forms a probability near 1
  log_F <- function(x) ifelse(x > 0, -log1p(exp(-x)), x - log1p(exp(x)))
  x <- c(-1000, -800, -40, -1, 0, 1, 40, 800, 1000)
  expect_equal(brms:::log_inv_logit(x), log_F(x))
  expect_equal(brms:::log1m_inv_logit(x), log_F(-x))
  expect_equal(brms:::log_inv_logit(-800), -800)
  expect_equal(brms:::log1m_inv_logit(40), -40 - log1p(exp(-40)))
  # the naive forms are accurate in between, and these agree with them there
  x <- c(-5, -0.5, 0, 0.5, 5)
  expect_equal(brms:::log_inv_logit(x), log(1 / (1 + exp(-x))))
  expect_equal(brms:::log1m_inv_logit(x), log(1 - 1 / (1 + exp(-x))))
  # Stan's conventions for non-finite input, and the shape of the input
  expect_true(is.nan(brms:::log_inv_logit(NaN)))
  expect_equal(brms:::log_inv_logit(c(-Inf, Inf)), c(-Inf, 0))
  m <- matrix(c(-50, 2, -0.5, 50), nrow = 2)
  expect_equal(dim(brms:::log1m_inv_logit(m)), c(2L, 2L))
})

test_that("log_softmax does not overflow", {
  # exp() overflows above 709.78, which made every category -Inf; #1401
  # routed softmax() through log_softmax() but did not shift it
  x <- matrix(c(0, 800, 801), nrow = 1)
  expect_equal(brms:::log_softmax(x), x - 801 - log1p(exp(-1)))
  expect_equal(sum(brms:::softmax(x)), 1)
  # draws x observations x categories, as inv_link_categorical() passes it
  a <- array(c(0, 900, 0, 1, 2, 3), dim = c(1, 2, 3))
  expect_equal(brms:::log_softmax(a)[1, 2, ], c(0, -899, -897))
  # it underflowed too, giving +Inf, when every category was below -745
  x <- matrix(c(-800, -801), nrow = 1)
  expect_equal(brms:::log_softmax(x), x + 800 - log1p(exp(-1)))
  # and it agrees with the unshifted form where that is accurate, including
  # a non-finite maximum, which is left unshifted
  x <- matrix(c(-1, 0.5, 2, 3, -2, 0), nrow = 2)
  expect_equal(brms:::log_softmax(x), x - log(rowSums(exp(x))))
  x <- matrix(c(Inf, 0), nrow = 1)
  expect_identical(brms:::log_softmax(x), x - log(rowSums(exp(x))))
  # the reshape keeps names, and zero categories pass through
  x <- array(1:12, dim = c(2, 2, 3))
  dimnames(x) <- list(NULL, NULL, c("a", "b", "c"))
  expect_identical(dimnames(brms:::log_softmax(x)), dimnames(x))
  expect_identical(dim(brms:::log_softmax(matrix(0, 4, 0))), c(4L, 0L))
})

test_that("the cloglog and softit log CDFs stay exact in the lower tail", {
  # both formed exp(x) first, which is subnormal below -708.4 and 0 below
  # -745.13; below -37 the log CDF is x to double precision
  x <- c(-800, -746, -744, -720, -40)
  expect_equal(brms:::log_cdf(x, "cloglog"), x)
  expect_equal(brms:::log_inv_softit(x), x)
  # above the cutoff the exact forms still apply
  expect_equal(brms:::log_cdf(-10, "cloglog"), log(-expm1(-exp(-10))))
  expect_equal(
    brms:::log_inv_softit(c(-10, 0, 3)),
    log(log1p(exp(c(-10, 0, 3)))) - log1p(log1p(exp(c(-10, 0, 3))))
  )
  # NaN stays NaN and empty input stays numeric
  expect_true(is.nan(brms:::log_cdf(NaN, "cloglog")))
  expect_identical(brms:::log_cdf(numeric(0), "cloglog"), numeric(0))
})
