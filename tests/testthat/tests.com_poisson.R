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

test_that("COM-Poisson log_lik applies censoring, truncation, and weights", {
  for (nu in c(1, 2)) {
    mu <- if (nu == 1) 2 else 0.8
    ref <- cmp_reference(mu, nu)
    for (bounds in list(c(-Inf, Inf), c(0, Inf), c(-Inf, 5), c(0, 5), c(0.2, 5))) {
      for (cens in c(0, -1, 1, 2)) for (weight in c(1, 2.5)) {
        prep <- structure(list(
          family = brmsfamily("com_poisson"),
          dpars = list(mu = c(mu, mu), shape = nu),
          data = list(Y = 1, cens = cens, rcens = 3,
                      weights = weight)), class = "brmsprep")
        if (is.finite(bounds[1])) prep$data$lb <- bounds[1]
        if (is.finite(bounds[2])) prep$data$ub <- bounds[2]
        mass <- switch(as.character(cens),
          "0" = ref$p[2], "-1" = sum(ref$p[ref$y <= 1]),
          "1" = sum(ref$p[ref$y > 1]),
          "2" = sum(ref$p[ref$y > 1 & ref$y <= 3]))
        norm <- sum(ref$p[ref$y >= ceiling(bounds[1]) & ref$y <= bounds[2]])
        expected <- weight * (log(mass) - log(norm))
        expect_equal(brms:::log_lik_com_poisson(1, prep), rep(expected, 2),
                     tolerance = 1e-10)
      }
    }
  }
})

test_that("COM-Poisson censored and truncated likelihoods match Stan", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
  chunk <- system.file("chunks", "fun_com_poisson.stan", package = "brms")
  code <- c("functions {", readLines(chunk), "}",
    "data { real mu; real nu; }",
    "generated quantities { vector[4] ll;",
    "real lnorm = log_diff_exp(com_poisson_lcdf(5 | mu, nu), com_poisson_lcdf(-1 | mu, nu));",
    "ll[1] = com_poisson_lpmf(1 | mu, nu) - lnorm;",
    "ll[2] = com_poisson_lcdf(1 | mu, nu) - lnorm;",
    "ll[3] = com_poisson_lccdf(1 | mu, nu) - lnorm;",
    "ll[4] = log_diff_exp(com_poisson_lcdf(3 | mu, nu), com_poisson_lcdf(1 | mu, nu)) - lnorm; }")
  sf <- tempfile(fileext = ".stan")
  writeLines(code, sf)
  mod <- cmdstanr::cmdstan_model(sf, quiet = TRUE)
  fit <- mod$sample(data = list(mu = 0.8, nu = 2), fixed_param = TRUE,
                   chains = 1, iter_sampling = 1, refresh = 0,
                   show_messages = FALSE, sig_figs = 18)
  actual <- as.numeric(posterior::as_draws_matrix(fit$draws("ll"))[1, ])
  prep <- structure(list(family = brmsfamily("com_poisson"),
                        dpars = list(mu = 0.8, shape = 2),
                        data = list(Y = 1, lb = 0, ub = 5, rcens = 3)),
                    class = "brmsprep")
  expected <- vapply(c(0, -1, 1, 2), function(cens) {
    prep$data$cens <- cens
    brms:::log_lik_com_poisson(1, prep)
  }, numeric(1))
  expect_equal(actual, expected, tolerance = 1e-8)
  sc <- make_stancode(y | trunc(lb = 0, ub = 5) ~ 1,
                     data = data.frame(y = 0:5), family = brmsfamily("com_poisson"))
  expect_match(sc, "com_poisson_lcdf(lb[n] - 1", fixed = TRUE)
})
