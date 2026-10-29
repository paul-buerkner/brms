context("Native Cox predictions")

cox_test_prep <- function(ns = 20L, nobs = 3L, tail = "constant") {
  basis <- brms:::bhaz_basis_matrix(
    seq(0, 4, length.out = 30),
    list(df = 5, intercept = TRUE, Boundary.knots = c(0, 4))
  )
  sb <- matrix(rep(c(0.1, 0.15, 0.2, 0.25, 0.3), each = ns), ns)
  fam <- cox()
  fam$fun <- "cox"
  structure(list(
    ndraws = ns, nobs = nobs, family = fam,
    dpars = list(mu = matrix(seq(0.5, 3, length.out = ns), ns, nobs)),
    data = list(Y = as.array(rep(1, nobs))),
    bhaz = list(
      basis = basis, sbhaz = list(sb), Jgrbhaz = rep(1L, nobs),
      version = package_version("2.23.2"), tail = tail
    )
  ), class = "brmsprep")
}

test_that("Cox Stan data use response-specific names", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, y2 = 2:11, x = seq_len(10))
  form <- bf(y ~ x, family = cox()) +
    bf(y2 ~ x, family = cox()) + set_rescor(FALSE)
  sdata <- standata(form, dat)
  expect_false(anyDuplicated(names(sdata)) > 0L)
  expect_true(all(c("Zbhaz_y", "Zcbhaz_y", "Kbhaz_y", "con_sbhaz_y",
                    "Zbhaz_y2", "Zcbhaz_y2") %in% names(sdata)))
  fit <- brm(form, dat, empty = TRUE)
  frame <- brms:::brmsframe(fit$formula, data = fit$data, basis = fit$basis)
  pars <- c("b_y_Intercept", "b_y_x", "b_y2_Intercept", "b_y2_x",
             paste0("sbhaz_y[", 1:5, "]"), paste0("sbhaz_y2[", 1:5, "]"))
  draws <- posterior::as_draws_matrix(
    matrix(0.2, 3, length(pars), dimnames = list(NULL, pars))
  )
  prep <- prepare_predictions(
    frame, draws = draws, sdata = sdata,
    cox_version = fit$version$brms, tail = "constant"
  )
  expect_equal(prep$resps$y$bhaz$basis, fit$basis$resps$y$bhaz$basis_matrix)
  expect_equal(prep$resps$y2$bhaz$basis, fit$basis$resps$y2$bhaz$basis_matrix)
  expect_equal(dim(posterior_predict(prep)), c(3L, 10L, 2L))
  expect_equal(dim(posterior_epred(
    prep, horizon = 2, dpar = NULL, nlpar = NULL, sort = FALSE
  )), c(3L, 10L, 2L))
})

test_that("Cox preparation preserves singleton and absent strata", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, x = seq_len(10),
                    g = factor(rep(c("a", "b"), 5)))
  fit <- brm(y | bhaz(gr = g) ~ x, data = dat, family = cox(), empty = TRUE)
  frame <- brms:::brmsframe(fit$formula, data = fit$data, basis = fit$basis)
  pars <- c(paste0("sbhaz[a,", 1:5, "]"), paste0("sbhaz[b,", 1:5, "]"))
  draws <- posterior::as_draws_matrix(
    matrix(0.2, 3, 10, dimnames = list(NULL, pars))
  )
  for (rows in list(1:2, 1L)) {
    sdata <- standata(fit, newdata = dat[rows, ])
    prep <- brms:::prepare_predictions_bhaz(
      frame$dpars$mu, draws, sdata, cox_version = fit$version$brms
    )
    expect_equal(dim(prep$bhaz), c(3L, length(rows)))
    expect_equal(prep$basis, fit$basis$bhaz$basis_matrix)
    expect_equal(prep$bhaz, unname(tcrossprod(matrix(0.2, 3, 5), sdata$Zbhaz)))
    expect_equal(prep$cbhaz, unname(tcrossprod(matrix(0.2, 3, 5), sdata$Zcbhaz)))
  }
})

test_that("Unchecked Cox responses do not extrapolate the spline", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, x = seq_len(10))
  fit <- brm(y ~ x, dat, family = cox(), empty = TRUE)
  expect_silent(sdata <- standata(
    fit, newdata = data.frame(y = c(0, 100), x = 1:2),
    check_response = FALSE
  ))
  expect_equal(as.vector(sdata$Zcbhaz %*% rep(0.2, 5)), c(0, 1))
  expect_silent(sdata <- standata(
    fit, newdata = data.frame(x = 1:2), check_response = FALSE
  ))
  expect_true(all(is.na(sdata$Y)))
  expect_equal(as.vector(sdata$Zcbhaz %*% rep(0.2, 5)), c(0, 0))
})

test_that("Cox preparation uses native links and nonlinear parameter dispatch", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, x = seq_len(10))
  for (link in c("log", "identity", "softplus", "squareplus")) {
    fit <- brm(y ~ 0 + Intercept, dat, family = cox(link), empty = TRUE)
    frame <- brms:::brmsframe(fit$formula, data = fit$data, basis = fit$basis)
    pars <- c("b_Intercept", paste0("sbhaz[", 1:5, "]"))
    draws <- posterior::as_draws_matrix(
      matrix(rep(c(2, rep(0.2, 5)), each = 3), 3,
             dimnames = list(NULL, pars))
    )
    prep <- prepare_predictions(
      frame, draws = draws, sdata = standata(fit),
      cox_version = fit$version$brms, tail = "constant"
    )
    expect_equal(brms:::get_dpar(prep, "mu"),
                 matrix(brms:::inv_link(2, link), 3, 10))
    expect_true(all(is.finite(posterior_predict(prep))))
  }
  fit <- brm(bf(y ~ exp(a), a ~ 1, nl = TRUE), dat,
             family = cox(), empty = TRUE)
  frame <- brms:::brmsframe(fit$formula, data = fit$data, basis = fit$basis)
  pars <- c("b_a_Intercept", paste0("sbhaz[", 1:5, "]"))
  draws <- posterior::as_draws_matrix(
    matrix(rep(c(0.2, rep(0.2, 5)), each = 3), 3,
           dimnames = list(NULL, pars))
  )
  prep <- prepare_predictions(
    frame, draws = draws, sdata = standata(fit),
    cox_version = fit$version$brms, tail = "constant"
  )
  expect_equal(prep$bhaz$basis, fit$basis$bhaz$basis_matrix)
  expect_true(all(is.finite(posterior_predict(prep))))
})

test_that("Cox inversion and tails agree with independent forward evaluation", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep()
  sb <- prep$bhaz$sbhaz[[1]]
  basis <- prep$bhaz$basis
  mu <- prep$dpars$mu[, 1]
  probs <- seq(0.01, 0.99, length.out = prep$ndraws)
  for (tail in c("constant", "loglinear")) {
    q <- brms:::qcox(probs, mu, sb, basis, tail)
    expect_equal(brms:::pcox(q, mu, sbhaz = sb, basis = basis, tail = tail),
                 probs, tolerance = 1e-8)
    at <- c(0.2, 1, 2, 3.9)
    for (t in at) {
      expected <- mu * as.vector(sb %*%
        as.vector(brms:::bhaz_basis_matrix(t, basis = basis, integrate = TRUE)))
      expect_equal(brms:::pcox(t, mu, sbhaz = sb, basis = basis,
                               log.p = TRUE, lower.tail = FALSE), -expected)
    }
    f <- function(x) brms:::dcox(
      x, mu = 1, sbhaz = sb[1, , drop = FALSE], basis = basis, tail = tail
    )
    expect_equal(integrate(f, 0, Inf)$value, 1, tolerance = 1e-6)
  }
  q <- brms:::qcox(0.99, mu, sb, basis, "none")
  expect_true(all(q < 0))
  expect_equal(cox_decode(q)$time, rep(4, length(mu)))
  expect_true(all(is.na(brms:::pcox(
    5, mu, sbhaz = sb, basis = basis, tail = "none"
  ))))
  expect_equal(brms:::qcox(0, mu, sb, basis), rep(0, length(mu)))
  expect_true(all(is.infinite(brms:::qcox(1, mu, sb, basis))))
  expect_warning(q <- brms:::qcox(c(-0.1, NA, rep(0.5, 18)), mu, sb, basis),
                 "NaNs produced")
  expect_true(all(is.na(q[1:2])))
  expect_equal(brms:::dcox(Inf, mu, sbhaz = sb, basis = basis,
                           tail = "loglinear"), rep(0, length(mu)))
  rate <- 1 / as.vector(brms:::bhaz_basis_matrix(4, basis = basis))[5] / 4
  flat <- matrix(c(rep((1 - rate) / 4, 4), rate), 1)
  expect_equal(brms:::hcox(Inf, 1, sbhaz = flat, basis = basis,
                           tail = "loglinear"), 0.25)
})

test_that("Native Cox outputs and stable truncation use the same distribution", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep()
  pred <- posterior_predict(prep)
  expect_equal(dim(pred), c(20L, 3L))
  expect_true(all(pred >= 0))
  for (output in c("probability", "pit", "density", "quantile")) {
    args <- list(prep, output = output, sort = FALSE, p = 0.5)
    got <- do.call(posterior_predict, args)
    expect_equal(dim(got), c(20L, 3L))
    expect_true(all(is.finite(got)))
  }
  prep$dpars$mu[] <- 100
  prep$data$lb <- rep(1, 3)
  prep$data$ub <- rep(2, 3)
  q <- posterior_predict(prep, output = "quantile", p = 0.5)
  expect_true(all(q > 1 & q < 2))
  expect_equal(unname(posterior_predict(prep, output = "probability", q = q[, 1])),
               matrix(0.5, 20, 3), tolerance = 1e-7)
  expect_true(all(is.finite(posterior_predict(
    prep, output = "density", q = q[, 1], log = TRUE
  ))))
  expect_equal(unname(posterior_predict(prep, output = "quantile", p = 1)),
               matrix(2, 20, 3), tolerance = 1e-8)
  prep$bhaz$tail <- "none"
  prep$data$ub[] <- 5
  expect_error(posterior_predict(prep), "truncation bounds")
  prep$data <- list(Y = rep(1, 3))
  prep$bhaz$version <- package_version("2.23.1")
  expect_error(posterior_predict(prep), "Please refit")
  prep$bhaz$version <- package_version("2.23.2")
  prep$bhaz$basis <- NULL
  expect_error(posterior_predict(prep), "Please refit")
  prep <- cox_test_prep()
  prep$family <- mixture(cox(), cox())
  prep$dpars <- list(mu1 = prep$dpars$mu, mu2 = prep$dpars$mu)
  prep$bhaz <- list(mu1 = prep$bhaz, mu2 = prep$bhaz)
  child <- brms:::pseudo_prep_for_mixture(prep, 1)
  expect_error(posterior_predict(child), "non-mixture")
})

test_that("Cox inversion controls probability error and narrow truncation", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep(ns = 3L, nobs = 1L)
  prep$dpars$mu[, 1] <- c(1e2, 1e6, 1e10)
  prep$bhaz$sbhaz[[1]][] <- rep(
    c(1e-10, 1e-6, 0.01, 0.1, 0.8899989999), each = 3
  )
  q <- posterior_predict(prep, output = "quantile", p = 0.5)
  expect_equal(as.vector(posterior_predict(
    prep, output = "probability", q = q[, 1]
  )), rep(0.5, 3), tolerance = 1e-8)

  prep <- cox_test_prep(nobs = 1L)
  prep$dpars$mu[] <- 1
  prep$data$lb <- 1.0259159159159161
  prep$data$ub <- prep$data$lb + 1e-10
  q <- posterior_predict(prep, output = "quantile", p = 0.5)
  expect_true(all(q >= prep$data$lb & q <= prep$data$ub))
  expect_equal(as.vector(posterior_predict(
    prep, output = "probability", q = q[, 1]
  )), rep(0.5, 20), tolerance = 1e-5)
  draws <- posterior_predict(prep)
  expect_true(all(draws >= prep$data$lb & draws <= prep$data$ub))
})

test_that("Cox quantiles retain extreme log and upper-tail probabilities", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep(nobs = 1L)
  prep$dpars$mu[] <- 1
  rate <- cox_bhaz(4, prep)[, 1]
  cases <- list(
    list(p = -100, lower.tail = FALSE, log.p = TRUE),
    list(p = -1e-20, lower.tail = TRUE, log.p = TRUE),
    list(p = 1e-20, lower.tail = FALSE, log.p = FALSE)
  )
  for (args in cases) {
    q <- do.call(posterior_predict,
                 c(list(prep, output = "quantile"), args))
    H <- do.call(qexp, args)
    expect_equal(as.vector(q), 4 + (H - 1) / rate, tolerance = 1e-10)
  }
  q <- posterior_predict(prep, output = "quantile", p = -500, log.p = TRUE)
  expect_true(all(is.finite(q) & q > 0))
  expect_equal(as.vector(posterior_predict(
    prep, output = "probability", q = q[, 1], log.p = TRUE
  )), rep(-500, 20), tolerance = 1e-8)
  prep$dpars$mu[] <- 1e4
  prep$data$lb <- 1
  prep$data$ub <- 3
  q <- posterior_predict(prep, output = "quantile", p = -1000,
                         lower.tail = FALSE, log.p = TRUE)
  expect_true(all(is.finite(q) & q >= 1 & q <= 3))
  expect_equal(as.vector(posterior_predict(
    prep, output = "probability", q = q[, 1], lower.tail = FALSE, log.p = TRUE
  )), rep(-1000, 20), tolerance = 1e-7)
})

test_that("Cox survival curves and streaming summaries agree", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep()
  times <- c(0, 1, 4)
  curves <- cox_surv_curve(times, prep)
  expect_equal(dim(curves), c(20L, 3L, 3L))
  expect_equal(curves[, 3, 1], exp(-prep$dpars$mu[, 1]))
  expect_equal(cox_surv(1, prep)[, 1], curves[, 2, 1])
  summary <- cox_surv_summary(times, prep, probs = c(0.1, 0.5))
  expect_equal(summary[, 1, 1], apply(curves[, , 1], 2, quantile, probs = 0.1))
  expect_equal(as.vector(cox_cumbhaz(c(0, 4), prep)[, 2]), rep(1, 20))
  prep$old_order <- c(3L, 1L, 2L)
  expect_equal(cox_surv(c(0.5, 1.5, 2.5), prep),
               exp(-prep$dpars$mu * cox_cumbhaz(c(0.5, 1.5, 2.5), prep)))
  expect_equal(cox_surv(1, prep), brms:::reorder_obs(cox_surv(1, prep,
    sort = TRUE), prep$old_order))
  expect_error(cox_surv_curve(numeric(0), prep), "times")
  expect_error(cox_surv(1:2, prep), "length")
})

test_that("Cox RMST agrees with independent integration and truncation", {
  skip_if_not_installed("splines2")
  skip_if_not_installed("statmod")
  prep <- cox_test_prep(ns = 3L, nobs = 1L)
  sb <- prep$bhaz$sbhaz[[1]][1, , drop = FALSE]
  basis <- prep$bhaz$basis
  for (horizon in c(0, 2, 6)) {
    got <- posterior_epred(prep, horizon = horizon, dpar = NULL,
                            nlpar = NULL, sort = FALSE)
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
  expect_equal(as.vector(posterior_epred(
    prep, horizon = 0.5, dpar = NULL, nlpar = NULL, sort = FALSE
  )), rep(0.5, 3))
  set.seed(42)
  pred <- replicate(2000, posterior_predict(prep)[1, 1])
  got <- posterior_epred(prep, horizon = 1.5, dpar = NULL,
                          nlpar = NULL, sort = FALSE)
  expect_equal(got[1, 1], mean(pmin(pred, 1.5)), tolerance = 0.02)
  expect_error(posterior_epred(prep, dpar = NULL, nlpar = NULL,
                               sort = FALSE), "horizon")
  prep$bhaz$tail <- "none"
  expect_error(posterior_epred(prep, horizon = 5, dpar = NULL,
                               nlpar = NULL, sort = FALSE), "support")
  prep <- cox_test_prep(ns = 2L)
  prep$data$lb <- c(0.25, 1, 2)
  prep$data$ub <- c(2.5, 3.5, Inf)
  prep$bhaz$tail <- "loglinear"
  got <- posterior_epred(prep, horizon = 6, dpar = NULL, nlpar = NULL,
                          sort = FALSE)
  for (i in seq_len(prep$nobs)) {
    ref <- integrate(function(t) brms:::pcox(
      t, prep$dpars$mu[1, i], sbhaz = prep$bhaz$sbhaz[[1]][1, , drop = FALSE],
      basis = prep$bhaz$basis, tail = "loglinear", lower.tail = FALSE,
      lb = prep$data$lb[i], ub = prep$data$ub[i]
    ), 0, 6, subdivisions = 1000, rel.tol = 1e-9)$value
    expect_equal(got[1, i], ref, tolerance = 1e-7)
  }
  prep$dpars$mu[] <- -1
  expect_error(posterior_epred(prep, horizon = 6, dpar = NULL, nlpar = NULL,
                               sort = FALSE), "finite positive")
})

test_that("Censored PIT randomizes after aggregation with shared uniforms", {
  skip_if_not_installed("splines2")
  prep <- cox_test_prep(nobs = 4L)
  prep$data$cens <- c(0, 1, -1, 2)
  prep$data$rcens <- rep(2, 4)
  F <- posterior_predict(prep, output = "probability")
  set.seed(41)
  avg <- cox_pit(prep)
  set.seed(41)
  conditional <- cox_pit(prep, average = FALSE)
  expect_equal(avg, unname(colMeans(conditional)))
  expect_equal(avg[1], mean(F[, 1]))
  expect_true(avg[2] >= mean(F[, 2]) && avg[2] <= 1)
  expect_true(avg[3] >= 0 && avg[3] <= mean(F[, 3]))
  upper <- mean(brms:::posterior_predict_cox(4, prep, output = "probability",
                                            q = 2))
  expect_true(avg[4] >= mean(F[, 4]) && avg[4] <= upper)
  prep$data$rcens <- NULL
  expect_error(cox_pit(prep), "rcens")
  prep$data$cens <- NULL
  prep$data$lb <- rep(0.5, 4)
  expect_equal(cox_snell(prep)[, 1],
    cox_cumhaz(1, prep)[, 1] - cox_cumhaz(0.5, prep)[, 1])
  prep$data$ub <- rep(2, 4)
  expect_error(cox_snell(prep), "upper truncation")
})

test_that("Censoring preserves existing signs and validates shapes", {
  draws <- matrix(c(1, -2, 3, 4), 2)
  expect_equal(cox_censor(draws, c(1.5, 3.5)),
               matrix(c(1, -1.5, 3, -3.5), 2))
  expect_equal(cox_censor(draws, Inf), draws)
  expect_error(cox_censor(draws, 1:3), "length")
  expect_error(cox_censor(draws, matrix(1, 1, 4)), "dimensions")
  expect_error(cox_censor(draws, 0), "positive")
  expect_equal(cox_censor_times(1:3, c(1, 1, 1), 5), rep(Inf, 5))
  skip_if_not_installed("prodlim")
  set.seed(12)
  ct <- cox_censor_times(1:4, c(1, 0, 0, 1), 1000)
  expect_true(all(ct %in% c(2, 3, Inf)))
  expect_true(all(c(2, 3, Inf) %in% ct))
})

test_that("Cox LOO-PIT rejects precomputed weights for selected or reordered draws", {
  skip_if_not_installed("splines2")
  testthat::local_mocked_bindings(
    prepare_cox_predictions = function(...) cox_test_prep(), .package = "brms"
  )
  fit <- brms:::brmsfit_example1
  ids <- seq_len(posterior::ndraws(fit))
  expect_error(loo_pit_cox(fit, loo_obj = list(), draw_ids = ids[-1]),
               "Recompute PSIS")
  expect_error(loo_pit_cox(fit, loo_obj = list(), draw_ids = rev(ids)),
               "reordered draws")
  expect_error(loo_pit_cox(fit, loo_obj = list(), newdata = fit$data),
               "newdata")
  wrong <- structure(list(psis_object = TRUE), yhash = "another-response")
  expect_error(loo_pit_cox(fit, loo_obj = wrong), "response")
})

test_that("Cox LOO-PIT computes likelihoods and CDFs from one preparation", {
  skip_if_not_installed("splines2")
  fit <- brms:::brmsfit_example1
  ns <- min(200L, posterior::ndraws(fit))
  prep <- cox_test_prep(ns = ns)
  set.seed(316)
  ll <- matrix(rnorm(ns * 3, -2, 0.1), ns)
  testthat::local_mocked_bindings(
    prepare_cox_predictions = function(...) prep,
    log_lik = function(object, ...) {
      stopifnot(identical(object, prep))
      ll
    },
    r_eff_log_lik = function(...) rep(1, 3),
    .package = "brms"
  )
  # The example fit has too few draws for reliable Pareto-k diagnostics.
  psis <- suppressWarnings(loo::psis(-ll, r_eff = rep(1, 3), cores = 1))
  expected <- colSums(
    weights(psis, normalize = TRUE, log = FALSE) *
      posterior_predict(prep, output = "probability")
  )
  got <- suppressWarnings(loo_pit_cox(
    fit, draw_ids = seq_len(ns), randomize_cens = FALSE,
    newdata = data.frame(g = "new"), allow_new_levels = TRUE, cores = 1
  ))
  expect_equal(got, expected)
})
