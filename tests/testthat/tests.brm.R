# calling context() avoids a strange bug in testthat 2.0.0
# cannot actually run brms models in tests as it takes way too long
context("Tests for brms error messages")

test_that("brm works fully with mock backend", {
  skip_on_cran()
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))

  # Positive control - forced error gets thrown and propagated
  expect_error(brm(y ~ x + (1|g), dat, backend = "mock",
                   stan_model_args = list(compile_error = "Test error")),
               "Test error")

  # Positive control - bad Stan code from stanvars gets an error
  # test passes but prints output that is somehow impossible to suppress
  # expect_error(suppressMessages(
  #    brm(y ~ x + (1|g), dat, backend = "mock",
  #        stanvars = stanvar(scode = "invalid;", block = "model"))
  # ))

  # Testing some models
  mock_fit <- brm(y ~ x + (1|g), dat, mock_fit = 1, backend = "mock", rename = FALSE)
  expect_equal(mock_fit$fit, 1)
})

test_that("brm(file = xx) works fully with mock backend", {
  skip_on_cran()
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))

  file <- tempfile(fileext = ".rds")
  mock_fit1 <- brm(y ~ x + (1|g), dat, mock_fit = "stored", backend = "mock",
                   rename = FALSE, file = file)
  expect_true(file.exists(file))

  mock_fit2 <- brm(y ~ x + (1|g), dat, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file)
  expect_equal(mock_fit2$fit, "stored")

  # In default settings, even using different data/model should result in the
  # model being loaded from file
  changed_data <- dat[1:8, ]
  mock_fit2 <- brm(y ~ x + 0, changed_data, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file)
  expect_equal(mock_fit2$fit, "stored")

  # Now test using file_refit = "on_change" which should be more clever
  # No change
  mock_fit2 <- brm(y ~ x + (1|g), dat, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file)
  expect_equal(mock_fit2$fit, "stored")


  # Change data, but not code
  mock_fit2 <- brm(y ~ x + (1|g), changed_data, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file, file_refit = "on_change")
  expect_equal(mock_fit2$fit, "new")

  # Change code but not data
  mock_fit2 <- brm(y ~ x + (1|g), dat, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file, file_refit = "on_change",
                   prior = prior(normal(0,2), class = sd))
  expect_equal(mock_fit2$fit, "new")

  # Change both
  mock_fit2 <- brm(y ~ x + 0, changed_data, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file, file_refit = "on_change")
  expect_equal(mock_fit2$fit, "new")
})


test_that("brm produces expected errors", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))

  # formula parsing
  expect_error(brm(~ x + (1|g), dat, file = "test"),
               "Response variable is missing")
  expect_error(brm(bf(y ~ a, nl = TRUE)),
               "No non-linear parameters specified")
  expect_error(brm(bf(y | se(sei) ~ x, sigma ~ x), dat),
               "Cannot predict or fix 'sigma' in this model")
  expect_error(brm(y | se(sei) ~ x, dat, family = weibull()),
               "Argument 'se' is not supported for family")
  expect_error(brm(y | se(sei) + se(sei2) ~ x, dat, family = gaussian()),
               "Each addition argument may only be defined once")
  expect_error(brm(y | abc(sei) ~ x, family = gaussian()),
               "The following addition terms are invalid:\n'abc(sei)'",
               fixed = TRUE)
  expect_error(brm(y | disp(sei) ~ x, dat, family = gaussian()),
               "The following addition terms are invalid:")
  expect_error(brm(bf(y ~ x, shape ~ x), family = gaussian()),
               "The parameter 'shape' is not a valid distributional")
  expect_error(brm(y ~ x + (1|abc|g/x), dat),
               "Can only combine group-level terms")
  expect_error(brm(y ~ x + (1|g) + (x|g), dat),
               "Duplicated group-level effects are not allowed")
  expect_error(brm(y~mo(g)*t2(x), dat), fixed = TRUE,
               "The term 'mo(g):t2(x)' is invalid")
  expect_error(brm(y~x*cs(g), dat), fixed = TRUE,
               "The term 'x:cs(g)' is invalid")
  expect_error(brm(y~me(x, 2 * g)*me(x, g), dat),
               "Variable 'x' is used in different calls to 'me'")
  expect_error(brm(y ~ 1 + set_rescor(TRUE), data = dat),
               "Function 'set_rescor' should not be part")

  # autocorrelation
  expect_error(brm(y ~ ar(x+y, g), dat),
               "Cannot coerce 'x \\+ y' to a single variable name")
  expect_error(brm(y ~ ar(gr = g1/g2), dat),
               "Illegal grouping term 'g1/g2'")
  expect_error(brm(y ~ ma(x), dat, poisson()),
               "Please set cov = TRUE")
  expect_error(brm(bf(y ~ 1) + arma(x), dat),
               "Autocorrelation terms can only be specified")

  # ordinal models
  expect_error(brm(rating ~ treat + (cs(period)|subject),
                   data = inhaler, family = categorical()),
               "Category specific effects are not supported")

  # families and links
  expect_error(brm(y ~ x, dat, family = poisson("inverse")),
               "'inverse' is not a supported link for family 'poisson'")
  expect_error(brm(y ~ x, dat, family = c("weibull", "sqrt")),
               "'sqrt' is not a supported link for family 'weibull'")
  expect_error(brm(y ~ x, dat, family = c("categorical", "probit")),
               "'probit' is not a supported link for family 'categorical'")
  expect_error(brm(y ~ x, dat, family = "ordinal"),
              "ordinal is not a supported family")
})

skip_if_no_stanli <- function() {
  skip_on_cran()
  skip_if_not_installed("stanli", minimum_version = stanli_min_version())
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
}

test_that("brm fits and post-processes a model with the stanli backend", {
  skip_if_no_stanli()
  set.seed(1234)
  dat <- data.frame(
    y = rnorm(40), x = rnorm(40), g = rep(1:5, each = 8)
  )
  fit <- brm(y ~ x + (1 | g), dat, backend = "stanli", chains = 2,
             iter = 600, seed = 1, silent = 2, refresh = 0)
  expect_s3_class(fit, "brmsfit")
  expect_s4_class(fit$fit, "stanfit")
  expect_equal(fit$backend, "stanli")
  expect_true(!is.null(fit$version$stanli))
  expect_equal(ndraws(fit), 600)
  expect_false(any(grepl("^z_1", variables(fit))))
  expect_true(all(c("b_Intercept", "b_x", "sd_g__Intercept", "sigma") %in%
                    variables(fit)))

  expect_equal(dim(posterior_predict(fit)), c(600, 40))
  expect_equal(dim(log_lik(fit)), c(600, 40))
  expect_s3_class(suppressWarnings(loo(fit)), "loo")
  expect_equal(nrow(elapsed_time(fit)), 2)
  expect_named(control_params(fit), c("adapt_delta", "max_treedepth", "metric",
                                     "adapt_engaged", "stepsize"))

  expect_false(needs_recompilation(fit))
  fit2 <- update(fit, newdata = dat[1:30, ], seed = 2, silent = 2, refresh = 0)
  expect_equal(nobs(fit2), 30)
  expect_equal(fit2$backend, "stanli")

  file <- tempfile(fileext = ".rds")
  saveRDS(fit, file)
  expect_equal(fixef(readRDS(file)), fixef(fit))
})

test_that("the stanli backend handles sampler arguments", {
  skip_if_no_stanli()
  set.seed(4321)
  dat <- data.frame(y = rnorm(30), x = rnorm(30))
  fit <- brm(y ~ x, dat, backend = "stanli", chains = 2, iter = 400,
             warmup = 150, thin = 2, init = 0, seed = 1, silent = 2,
             refresh = 0, control = list(adapt_delta = 0.9, max_treedepth = 8))
  expect_equal(ndraws(fit), 250)
  expect_equal(control_params(fit)$adapt_delta, 0.9)
  expect_equal(control_params(fit)$max_treedepth, 8)
  expect_warning(
    fit <- brm(y ~ x, dat, backend = "stanli", chains = 1, iter = 300,
               seed = 1, silent = 2, refresh = 0,
               control = list(metric = "dense_e", stepsize = 0.5,
                              kappa = 0.7)),
    "ignored by backend 'stanli': 'kappa'$"
  )
  expect_equal(control_params(fit)$metric, "dense_e")
  expect_equal(control_params(fit)$stepsize, 0.5)
  fit <- update(fit, newdata = dat[1:20, ], silent = 2, refresh = 0)
  expect_equal(control_params(fit)$metric, "dense_e")

  inits <- list(
    list(Intercept = 0, b = array(0.1), sigma = 1),
    list(sigma = 2)
  )
  fit <- brm(y ~ x, dat, backend = "stanli", chains = 2, iter = 300,
             init = inits, seed = 1, silent = 2, refresh = 0)
  expect_equal(ndraws(fit), 300)
  expect_error(
    brm(y ~ x, dat, backend = "stanli", chains = 1, iter = 300,
        init = list(list(sigma = -1)), silent = 2, refresh = 0),
    "sigma"
  )

  empty <- brm(y ~ x, dat, backend = "stanli", chains = 0, silent = 2)
  expect_error(contains_draws(empty))
  filled <- update(empty, chains = 2, iter = 300, seed = 1, silent = 2,
                   refresh = 0, recompile = FALSE)
  expect_equal(ndraws(filled), 300)
})

test_that("log_prob-based methods need no compiler with the stanli backend", {
  skip_if_no_stanli()
  skip_if_not_installed("bridgesampling")
  set.seed(2468)
  dat <- data.frame(y = rnorm(40), x = rnorm(40))
  fit <- brm(y ~ x, dat, backend = "stanli", chains = 2, iter = 1000,
             seed = 1, silent = 2, refresh = 0,
             save_pars = save_pars(all = TRUE))
  expect_no_message(
    ml <- bridge_sampler(fit, silent = TRUE), message = "Recompiling"
  )
  expect_true(is.finite(ml$logml))
  file <- tempfile(fileext = ".rds")
  saveRDS(fit, file)
  ml2 <- bridge_sampler(readRDS(file), silent = TRUE)
  expect_equal(ml2$logml, ml$logml, tolerance = 0.05)
  expect_s3_class(suppressWarnings(loo(fit, moment_match = TRUE)), "loo")

  fits <- brm_multiple(y ~ x, list(dat, dat[1:30, ]), backend = "stanli",
                       chains = 1, iter = 1000, seed = 1, silent = 2,
                       refresh = 0, save_pars = save_pars(all = TRUE))
  expect_true(is.finite(bridge_sampler(fits, silent = TRUE)$logml))
})

test_that("the stanli backend runs the approximate algorithms", {
  skip_if_no_stanli()
  set.seed(1357)
  dat <- data.frame(
    y = rnorm(60), x = rnorm(60), g = rep(1:6, each = 10)
  )
  dat$y <- 0.5 + 0.8 * dat$x + dat$y
  ref <- brm(y ~ x + (1 | g), dat, backend = "stanli", chains = 2,
             iter = 1000, seed = 1, silent = 2, refresh = 0)

  for (algorithm in c("meanfield", "fullrank", "pathfinder", "laplace")) {
    fit <- suppressMessages(suppressWarnings(brm(
      y ~ x + (1 | g), dat, backend = "stanli", algorithm = algorithm,
      seed = 1, silent = 2, refresh = 0
    )))
    expect_s3_class(fit, "brmsfit")
    expect_equal(fit$algorithm, algorithm, info = algorithm)
    expect_equal(fit$backend, "stanli")
    expect_equal(nchains(fit), 1, info = algorithm)
    expect_equal(ndraws(fit), 1000, info = algorithm)
    expect_true(all(c("lp__", "lp_approx__") %in% variables(fit)),
                info = algorithm)
    expect_false(any(grepl("^z_1", variables(fit))))

    smry <- suppressWarnings(summary(fit))
    method <- if (algorithm %in% c("meanfield", "fullrank")) {
      paste0("variational(", algorithm, ")")
    } else {
      paste0(algorithm, "(lbfgs)")
    }
    expect_equal(smry$sampler, method, info = algorithm)
    expect_equal(smry$total_ndraws, 1000)
    expect_output(print(smry), "Draws were sampled using", fixed = TRUE)
    expect_equal(dim(posterior_predict(fit)), c(1000, 60))
    expect_equal(dim(posterior_epred(fit)), c(1000, 60))
    expect_equal(dim(log_lik(fit)), c(1000, 60))
    draws <- as_draws_df(fit)
    expect_equal(nrow(draws), 1000)
    expect_true(all(c("b_Intercept", "b_x", "sd_g__Intercept", "sigma",
                      "lp__", "lp_approx__") %in% names(draws)))
    expect_s3_class(suppressWarnings(loo(fit)), "loo")
    expect_equal(nrow(elapsed_time(fit)), 1)

    expect_equal(
      fixef(fit)["x", "Estimate"], fixef(ref)["x", "Estimate"],
      tolerance = 3 * fixef(ref)["x", "Est.Error"], scale = 1,
      info = algorithm
    )
    expect_error(nuts_params(fit))
    expect_false(needs_recompilation(fit))
  }

  fit <- suppressWarnings(brm(
    y ~ x, dat, backend = "stanli", algorithm = "meanfield", iter = 3000,
    seed = 2, silent = 2, refresh = 0, draws = 500, init = 0,
    control = list(tol_rel_obj = 0.005)
  ))
  expect_equal(ndraws(fit), 500)
  expect_message(
    fit <- brm(y ~ x, dat, backend = "stanli", algorithm = "pathfinder",
               chains = 4, seed = 2, silent = 1, refresh = 0, draws = 300,
               init = list(list(Intercept = 0, sigma = 1))),
    "single-path Pathfinder only. Running 1 path instead of 4"
  )
  expect_equal(ndraws(fit), 300)
  expect_equal(nchains(fit), 1)
  expect_message(
    fit1 <- brm(y ~ x, dat, backend = "stanli", algorithm = "pathfinder",
                chains = 1, seed = 2, silent = 2, refresh = 0, draws = 300,
                init = list(list(Intercept = 0, sigma = 1))),
    NA
  )
  expect_equal(fixef(fit1), fixef(fit))
  fit <- brm(y ~ x, dat, backend = "stanli", algorithm = "laplace",
             seed = 2, silent = 2, refresh = 0, draws = 200,
             init = list(list(sigma = 1)))
  expect_equal(ndraws(fit), 200)
  fit2 <- update(fit, newdata = dat[1:40, ], seed = 3, silent = 2,
                 refresh = 0)
  expect_equal(fit2$algorithm, "laplace")
  expect_equal(nobs(fit2), 40)
  expect_equal(ndraws(fit2), 200)
  expect_error(
    brm(y ~ x, dat, backend = "stanli", algorithm = "laplace", silent = 2,
        refresh = 0, init = list(list(sigma = 1), list(sigma = 2))),
    "Algorithm 'laplace' takes 1"
  )
})

test_that("the stanli backend runs fixed_param", {
  skip_if_no_stanli()
  set.seed(2468)
  dat <- data.frame(y = rnorm(30), x = rnorm(30))
  inits <- replicate(2, list(Intercept = 1, b = array(0.5), sigma = 2),
                     simplify = FALSE)
  fit <- brm(y ~ x, dat, backend = "stanli", algorithm = "fixed_param",
             chains = 2, iter = 300, warmup = 100, init = inits, seed = 1,
             silent = 2, refresh = 0)
  expect_equal(fit$algorithm, "fixed_param")
  expect_equal(ndraws(fit), 400)
  draws <- as_draws_df(fit)
  expect_true(all(draws$b_x == 0.5))
  expect_true(all(draws$sigma == 2))
  expect_equal(summary(fit)$sampler, "sampling(Fixed_param)")
  expect_equal(dim(posterior_predict(fit)), c(400, 30))
  expect_true(sd(posterior_predict(fit)[, 1]) > 1)
  expect_equal(nrow(elapsed_time(fit)), 2)
})

test_that("the stanli backend refuses what it does not support", {
  skip_if_no_stanli()
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_error(
    brm(y ~ x, dat, backend = "stanli", opencl = opencl(c(0, 0)), silent = 2),
    "OpenCL is not supported by backend 'stanli'"
  )
  expect_error(
    brm(y ~ x, dat, backend = "stanli", future = TRUE, silent = 2),
    "'future' is not supported by backend 'stanli'"
  )
  expect_error(
    brm(y ~ x, dat, backend = "stanli", silent = 2,
        stanvars = stanvar(scode = "real zz = ;", block = "tparameters")),
    "Syntax error"
  )
})

