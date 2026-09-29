context("COM-Poisson probability calculations")

# Independent reference, deliberately using a fixed and oversized support.
cmp_reference <- function(mu, shape, K = 2000L) {
  y <- 0:K
  lw <- shape * (y * log(mu) - lgamma(y + 1))
  w <- exp(lw - max(lw))
  p <- w / sum(w)
  list(y = y, p = p, cdf = cumsum(p), mean = sum(y * p),
       logZ = max(lw) + log(sum(w)))
}

test_that("COM-Poisson CDF vectorization respects elementwise parameters", {
  x <- c(0, 1, 2, 4)
  mu <- c(1, 2, 3, 4)
  expect_warning(p <- brms:::pcom_poisson(x, mu, 1), NA)
  expect_equal(p, ppois(x, mu), tolerance = 1e-12)
  for (x in list(c(0, 1, 2, 4), c(1, 2, 4, 3))) {
    mu <- c(0.4, 0.8, 1.1, 1.2)
    shape <- c(0.5, 1.5, 3, 0.8)
    ref <- mapply(function(y, m, s) sum(cmp_reference(m, s)$p[seq_len(y + 1)]),
                  x, mu, shape)
    expect_equal(brms:::pcom_poisson(x, mu, shape), ref, tolerance = 1e-12)
    expect_equal(brms:::pcom_poisson(x, mu, shape),
                 mapply(brms:::pcom_poisson, x, mu, shape))
    xm <- matrix(x, nrow = 2)
    mm <- matrix(mu, nrow = 2)
    expect_equal(brms:::pcom_poisson(xm, mm, shape), matrix(ref, 2),
                 tolerance = 1e-12)
  }
  expect_equal(brms:::pcom_poisson(0:4, 2, 1), ppois(0:4, 2))
})
