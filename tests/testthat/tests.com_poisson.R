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

test_that("COM-Poisson CDF handles support boundaries and tail flags", {
  x <- c(-Inf, -2, -0.1, 0, 0.9, 1.8, Inf, NA_real_, NaN)
  for (lower in c(TRUE, FALSE)) for (lp in c(TRUE, FALSE)) {
    expect_equal(brms:::pcom_poisson(x, 2, 1, lower, lp),
                 ppois(x, 2, lower.tail = lower, log.p = lp))
    expect_equal(brms:::pcom_poisson(matrix(x, 3), 2, 1, lower, lp),
                 matrix(ppois(x, 2, lower.tail = lower, log.p = lp), 3))
  }
  expect_equal(brms:::pcom_poisson(c(-1, Inf), 0.8, 2), c(0, 1))
  set.seed(1938)
  u <- runif(1)
  set.seed(1938)
  pit <- brms:::pp_cdf(0, "com_poisson", NULL, NULL, randomized = TRUE,
                       mu = 2, shape = 1, int_response = TRUE)
  expect_equal(pit, u * ppois(0, 2))
  expect_equal(brms:::pp_cdf(1, "com_poisson", 0, 3, randomized = FALSE,
                             mu = 2, shape = 1, int_response = TRUE),
               ppois(1, 2) / ppois(3, 2))
})

test_that("Stan COM-Poisson CDF accepts negative cutoffs", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
  chunk <- system.file("chunks", "fun_com_poisson.stan", package = "brms")
  code <- c("functions {", readLines(chunk), "}",
    "data { int N; array[N] int y; vector[N] mu; vector[N] nu; }",
    "generated quantities { vector[N] lcdf; vector[N] lccdf;",
    "for (n in 1:N) { lcdf[n] = com_poisson_lcdf(y[n] | mu[n], nu[n]);",
    "lccdf[n] = com_poisson_lccdf(y[n] | mu[n], nu[n]); } }")
  sf <- tempfile(fileext = ".stan")
  writeLines(code, sf)
  mod <- cmdstanr::cmdstan_model(sf, quiet = TRUE)
  d <- expand.grid(y = c(-2L, -1L, 0L, 1L, 4L), nu = c(0.7, 1, 2))
  fit <- mod$sample(data = list(N = nrow(d), y = d$y, mu = rep(0.8, nrow(d)),
                                nu = d$nu), fixed_param = TRUE, chains = 1,
                   iter_sampling = 1, refresh = 0, show_messages = FALSE,
                   sig_figs = 18)
  draws <- posterior::as_draws_matrix(fit$draws())
  for (lower in c(TRUE, FALSE)) {
    variable <- if (lower) "lcdf" else "lccdf"
    expect_equal(as.numeric(draws[1, sprintf("%s[%d]", variable, seq_len(nrow(d)))]),
                 brms:::pcom_poisson(d$y, 0.8, d$nu, lower.tail = lower, log.p = TRUE),
                 tolerance = 1e-8)
  }
})
