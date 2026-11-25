# Unit test for paul-buerkner/brms/tests/testthat/test-add_rstan_model.R
# Verifies that add_rstan_model() synchronizes @model_pars and @par_dims on
# cmdstanr-backend brmsfit objects with excluded Stan parameters (e.g., r_1, Cor_1),
# preventing silent parameter corruption in .update_pars() and loo_moment_match().

test_that("add_rstan_model syncs @model_pars and @par_dims for cmdstanr fits", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rstan")
  skip_if_not_installed("loo")
  skip_on_cran()

  set.seed(2026)
  N <- 25
  T_w <- 4
  dat <- expand.grid(time = seq_len(T_w), id = seq_len(N))
  dat$x <- rnorm(nrow(dat))
  u0 <- rnorm(N, 0, 1.0)
  u1 <- 0.5 * u0 + rnorm(N, 0, 0.8)
  dat$y <- 1 + 0.8 * dat$x + u0[dat$id] + u1[dat$id] * dat$x +
    rnorm(nrow(dat), 0, 0.5)
  dat$y[1] <- dat$y[1] + 10

  fit <- brm(
    y ~ x + (1 + x | id),
    data = dat,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE),
    chains = 2, iter = 600, warmup = 300, cores = 2,
    seed = 123, refresh = 0, silent = 2
  )

  # cmdstanr omits excluded intermediate arrays (r_1, Cor_1) from @model_pars
  x <- suppressMessages(brms:::add_rstan_model(fit, overwrite = TRUE))
  expect_true(all(c("r_1", "Cor_1") %in% x$fit@model_pars))
  expect_identical(names(x$fit@par_dims), x$fit@model_pars)

  # Identity round-trip through unconstrain_pars -> .update_pars must not shift
  # downstream parameters (lprior, b_Intercept, cor_id__Intercept__x)
  dr0 <- as.matrix(x)
  sub_idx <- seq_len(min(50L, nrow(dr0)))
  up <- brms:::.unconstrain_pars(x, pars = dr0[sub_idx, , drop = FALSE])
  dr1 <- as.matrix(brms:::.update_pars(x, upars = up))

  expect_equal(dr1[, "b_Intercept"], dr0[sub_idx, "b_Intercept"], tolerance = 1e-4)
  expect_equal(dr1[, "cor_id__Intercept__x"], dr0[sub_idx, "cor_id__Intercept__x"], tolerance = 1e-4)
  expect_equal(dr1[, "lprior"], dr0[sub_idx, "lprior"], tolerance = 1e-4)
  expect_true(all(abs(dr1[, "cor_id__Intercept__x"]) <= 1))

  # Standard post-processing methods remain unaffected after slot synchronization
  expect_s3_class(summary(x), "brmssummary")
  expect_equal(posterior::as_draws_df(x), posterior::as_draws_df(fit))
  expect_equal(posterior_epred(x), posterior_epred(fit))
  expect_equal(log_lik(x), log_lik(fit))

  # loo_moment_match runs cleanly and improves high-k observations
  l <- suppressWarnings(loo(x))
  expect_gt(sum(l$diagnostics$pareto_k > 0.7), 0L)
  mm <- suppressWarnings(loo_moment_match(x, loo = l))
  expect_s3_class(mm, "psis_loo")
  expect_lt(sum(mm$diagnostics$pareto_k > 0.7), sum(l$diagnostics$pareto_k > 0.7))
})
