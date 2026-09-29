context("COM-Poisson probability calculations")

# Independent reference using log weights on a fixed, oversized support.
# Certify that omitted mass and first moment are negligible for every case.
cmp_reference <- function(mu, shape, K = 2000L) {
  y <- 0:K
  lw <- shape * (y * log(mu) - lgamma(y + 1))
  w <- exp(lw - max(lw))
  p <- w / sum(w)
  logZ <- max(lw) + log(sum(w))
  mean <- sum(y * p)
  ratio <- exp(shape * (log(mu) - log(K + 2)))
  log_tail <- shape * ((K + 1) * log(mu) - lgamma(K + 2)) -
    log1p(-ratio) - logZ
  log_moment_tail <- log_tail + log(K + 1 + ratio / (1 - ratio))
  stopifnot(ratio < 1, log_tail < log(1e-14),
            log_moment_tail < log(1e-14) + log(mean))
  list(y = y, p = p, cdf = cumsum(p), mean = mean, logZ = logZ)
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

test_that("COM-Poisson normalization, moments and tails agree with independent sums", {
  grid <- expand.grid(mu = c(0.05, 0.5, 1, 1.5, 2, 3.4, 10, 100, 1000),
                      nu = c(0.05, 0.2, 0.6, 1, 1.5, 5, 30, 100))
  refs <- lapply(seq_len(nrow(grid)), function(i) {
    cmp_reference(grid$mu[i], grid$nu[i], K = 20000)
  })
  expect_equal(brms:::log_Z_com_poisson(log(grid$mu), grid$nu),
               vapply(refs, `[[`, numeric(1), "logZ"), tolerance = 1e-10)
  expected_mean <- vapply(refs, `[[`, numeric(1), "mean")
  actual_mean <- brms:::mean_com_poisson(grid$mu, grid$nu)
  expect_true(all(abs(actual_mean - expected_mean) <= 1e-10 * pmax(1, expected_mean)))
  for (mu in c(0.05, 1.5, 2, 10)) for (nu in c(0.2, 1, 5, 30, 100)) {
    ref <- cmp_reference(mu, nu)
    expect_equal(sum(brms:::dcom_poisson(0:2000, mu, nu)), 1, tolerance = 1e-10)
    x <- 0:20
    p <- brms:::pcom_poisson(x, mu, nu)
    expect_equal(diff(c(0, p)), brms:::dcom_poisson(x, mu, nu), tolerance = 1e-10)
    expect_equal(p, ref$cdf[x + 1], tolerance = 1e-10)
    probs <- c(1e-8, 0.01, 0.5, 0.75, 0.99, 1 - 1e-8)
    q <- brms:::qcom_poisson(probs, mu, nu)
    qref <- vapply(probs, function(p) which(ref$cdf >= p)[1] - 1, numeric(1))
    # At a CDF jump indistinguishable from p in double precision (e.g. the
    # two almost equal masses at mu=2, shape=100), use the defining
    # inequalities instead of treating a rounded reference as exact.
    resolved <- abs(ref$cdf[qref + 1] - probs) > 1e-12
    expect_equal(q[resolved], qref[resolved])
    expect_true(all(brms:::pcom_poisson(q, mu, nu) >= probs - 1e-12))
    expect_true(all(brms:::pcom_poisson(q - 1, mu, nu) < probs + 1e-12))
  }
  # Survival probability is much smaller than machine epsilon, but its log
  # must remain finite for right censoring and far-tail truncation.
  y <- 6:100
  lw <- 30 * (y * log(2) - lgamma(y + 1))
  tail_ref <- max(lw) + log(sum(exp(lw - max(lw)))) - cmp_reference(2, 30)$logZ
  expect_equal(brms:::pcom_poisson(5, 2, 30, FALSE, TRUE), tail_ref, tolerance = 1e-10)
  prep <- structure(list(family = brmsfamily("com_poisson"),
    dpars = list(mu = 2, shape = 30), data = list(Y = 6, lb = 6, ub = 8)),
    class = "brmsprep")
  lp <- 30 * ((6:8) * log(2) - lgamma(7:9))
  expected <- lp[1] - (max(lp) + log(sum(exp(lp - max(lp)))))
  expect_equal(brms:::log_lik_com_poisson(1, prep), expected, tolerance = 1e-10)
  set.seed(1938)
  draws <- brms:::rcom_poisson(10000, 2, 30)
  expect_equal(mean(draws == 1), 0.5, tolerance = 0.02)
  expect_equal(mean(draws <= 2), brms:::pcom_poisson(2, 2, 30), tolerance = 0.001)
  expect_error(brms:::log_Z_com_poisson(log(2), 30, M = 0), "failed to converge")
  expect_error(brms:::mean_com_poisson(2, 30, M = 0), "failed to converge")
  expect_warning(q <- brms:::qcom_poisson(0.9, 2, 30, M = 1), "failed")
  expect_true(is.na(q))
  expect_equal(brms:::qcom_poisson(0.75, 2, 30, M = 2), 2)
})

test_that("compiled COM-Poisson values and gradients match independent references", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
  skip_if(cmdstanr::cmdstan_version() < "2.36.0")
  chunk <- system.file("chunks", "fun_com_poisson.stan", package = "brms")
  code <- c("functions {", readLines(chunk), "}",
    "data { int N; array[N] int y; vector[N] mu; vector[N] nu; int cap; int kind; int cutoff; }",
    "parameters { vector[2] theta; }",
    "model { if (cap == 0) target += com_poisson_sum(theta[1], exp(theta[2]), cap, 1e-12)[1];",
    "else if (kind == 0) target += com_poisson_log_lpmf(cutoff | theta[1], exp(theta[2]));",
    "else if (kind == 1) target += com_poisson_lcdf(cutoff | exp(theta[1]), exp(theta[2]));",
    "else target += com_poisson_lccdf(cutoff | exp(theta[1]), exp(theta[2])); }",
    "generated quantities { vector[N] lpmf; vector[N] lcdf; vector[N] lccdf; vector[N] means;",
    "for (n in 1:N) { lpmf[n] = com_poisson_lpmf(y[n] | mu[n], nu[n]);",
    "lcdf[n] = com_poisson_lcdf(y[n] | mu[n], nu[n]);",
    "lccdf[n] = com_poisson_lccdf(y[n] | mu[n], nu[n]);",
    "means[n] = com_poisson_sum(log(mu[n]), nu[n], cap, 1e-12)[2]; } }")
  sf <- tempfile(fileext = ".stan")
  writeLines(code, sf)
  mod <- cmdstanr::cmdstan_model(sf, quiet = TRUE)
  g <- expand.grid(mu = c(0.05, 1.5, 2, 10, 1000), nu = c(0.05, 0.6, 1, 5, 30, 100))
  d <- g[rep(seq_len(nrow(g)), each = 3), ]
  d$y <- as.integer(pmax(0, floor(d$mu) + rep(c(-1, 0, 5), nrow(g))))
  data <- list(N = nrow(d), y = d$y, mu = d$mu, nu = d$nu, cap = 10000L,
               kind = 0L, cutoff = 2L)
  fit <- mod$sample(data = data, init = list(list(theta = c(log(2), 0))),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18)
  draws <- posterior::as_draws_matrix(fit$draws())
  for (what in c("lpmf", "lcdf", "lccdf", "means")) {
    ref <- vapply(seq_len(nrow(d)), function(i) {
      r <- cmp_reference(d$mu[i], d$nu[i], K = 20000)
      if (what == "means") return(r$mean)
      y <- r$y
      lw <- d$nu[i] * (y * log(d$mu[i]) - lgamma(y + 1))
      if (what == "lpmf") return(lw[d$y[i] + 1] - r$logZ)
      lw <- if (what == "lcdf") lw[y <= d$y[i]] else lw[y > d$y[i]]
      max(lw) + log(sum(exp(lw - max(lw)))) - r$logZ
    }, numeric(1))
    actual <- as.numeric(draws[1, sprintf("%s[%d]", what, seq_len(nrow(d)))])
    expect_true(all(abs(actual - ref) <= 1e-8), info = what)
  }
  # CmdStan log_prob evaluates autodiff without sampling, including exactly
  # at shape=1 and either side of modal/old approximation switch boundaries.
  grad_grid <- expand.grid(mu = c(0.5, 1 - 1e-8, 1, 1 + 1e-8,
                                  1.5 - 1e-8, 1.5, 1.5 + 1e-8,
                                  2 - 1e-8, 2, 2 + 1e-8, 100),
                           nu = c(0.2, 1 - 1e-8, 1, 1 + 1e-8, 30))
  uf <- tempfile(fileext = ".json")
  df <- tempfile(fileext = ".json")
  of <- tempfile(fileext = ".csv")
  # Conditional moment identities also check derivatives through the CDF
  # tail selection and stable complement, as used by censored likelihoods.
  # Use one parameter vector per invocation: older CmdStan versions reshape
  # multi-row JSON parameter inputs in a different order.
  for (kind in 0:2) {
    data$kind <- kind
    cmdstanr::write_stan_json(data, df)
    gradients <- lapply(seq_len(nrow(grad_grid)), function(i) {
      cmdstanr::write_stan_json(list(params_r = c(log(grad_grid$mu[i]),
                                                 log(grad_grid$nu[i]))), uf)
      result <- processx::run(mod$exe_file(), c("log_prob", paste0("unconstrained_params=", uf),
        "jacobian=false", "data", paste0("file=", df), "output", paste0("file=", of),
        "sig_figs=18"), error_on_status = FALSE)
      expect_equal(result$status, 0L, info = result$stderr)
      read.csv(of, comment.char = "#")
    })
    gradients <- do.call(rbind, gradients)
    expected <- t(vapply(seq_len(nrow(grad_grid)), function(i) {
      mu <- grad_grid$mu[i]; nu <- grad_grid$nu[i]
      r <- cmp_reference(mu, nu, K = 20000)
      score <- r$y * log(mu) - lgamma(r$y + 1)
      keep <- switch(kind + 1L, r$y == data$cutoff,
                     r$y <= data$cutoff, r$y > data$cutoff)
      lw <- nu * score[keep]
      conditional <- exp(lw - max(lw))
      conditional <- conditional / sum(conditional)
      c(nu * (sum(conditional * r$y[keep]) - r$mean),
        nu * (sum(conditional * score[keep]) - sum(r$p * score)))
    }, numeric(2)))
    expect_true(all(abs(as.matrix(gradients[, 2:3]) - expected) <= 1e-7),
                info = paste("gradient kind", kind))
  }
  # An intentionally inadequate term limit must fail, not return a partial sum.
  data$cap <- 0L
  cmdstanr::write_stan_json(data, df)
  result <- processx::run(mod$exe_file(), c("log_prob", paste0("unconstrained_params=", uf),
    "jacobian=false", "data", paste0("file=", df), "output", paste0("file=", of)),
    error_on_status = FALSE)
  expect_true(result$status != 0L)
  expect_match(paste(result$stdout, result$stderr), "failed to converge")
})
