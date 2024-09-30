source("setup_tests_local.R")

test_that("Poisson model from brm doc works correctly", suppressWarnings({
  ## Poisson regression for the number of seizures in epileptic patients
  ## using student_t priors for population-level effects
  ## and half cauchy priors for standard deviations of group-level effects
  fit1 <- brm(
    count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
    data = epilepsy, family = poisson(),
    prior = prior(student_t(5,0,10), class = b) +
      prior(cauchy(0,2), class = sd),
    save_pars = save_pars(all = TRUE), refresh = 0,
    backend = "rstan"
  )
  print(fit1)
  ## generate a summary of the results
  expect_range(fixef(fit1)[4, 1], -0.28, -0.23)
  ## extract random effects standard devations and covariance matrices
  expect_equal(names(VarCorr(fit1)), c("obs", "patient"))
  ## extract group specific effects of each level
  expect_equal(names(ranef(fit1)), c("obs", "patient"))
  ## predict responses based on the fitted model
  expect_equal(dim(predict(fit1)), c(nobs(fit1), 4))
  ## plot conditional effects of each predictor
  me1 <- conditional_effects(fit1)
  expect_ggplot(plot(me1, ask = FALSE)[[4]])
  ## investigate model fit
  expect_range(WAIC(fit1)$estimates[3, 1], 1120, 1160)
  expect_ggplot(pp_check(fit1))

  # test kfold
  kfold1 <- kfold(fit1, chains = 1, iter = 1000, save_fits = TRUE)
  expect_range(kfold1$estimates[3, 1], 1210, 1260)
  # define a loss function
  rmse <- function(y, yrep) {
    yrep_mean <- colMeans(yrep)
    sqrt(mean((yrep_mean - y)^2))
  }
  # predict responses and evaluate the loss
  kfp1 <- kfold_predict(kfold1)
  rmse1 <- rmse(y = kfp1$y, yrep = kfp1$yrep)
  expect_range(rmse1, 6, 7)

  # test loo_moment_match
  loo1 <- loo(fit1)
  mmloo1 <- loo_moment_match(fit1, loo1, k_threshold = 1, cores = 1)
  expect_is(mmloo1, "loo")
}))

test_that("Ordinal model from brm doc works correctly", suppressWarnings({
  ## Ordinal regression modeling patient's rating of inhaler instructions
  ## category specific effects are estimated for variable 'treat'
  fit2 <- brm(
    rating ~ period + carry + cs(treat),
    data = inhaler, family = sratio("cloglog"),
    prior = set_prior("normal(0,5)"),
    iter = 1000, chains = 2, refresh = 0
  )
  print(fit2)
  expect_range(WAIC(fit2)$estimates[3, 1], 900, 950)
  suppressWarnings(ce <- conditional_effects(fit2, effect = "treat"),
                 "Predictions are treated as continuous variables")
  expect_ggplot(plot(ce)[[1]])
}))

test_that("Survival model from brm doc works correctly", suppressWarnings({
  ## Survival regression modeling the time between the first
  ## and second recurrence of an infection in kidney patients.
  fit3 <- brm(
    time | cens(censored) ~ age * sex + disease + (1|patient),
    data = kidney, family = lognormal(), refresh = 0,
    threads = threading(2, grainsize = 100),
    backend = "cmdstanr", save_pars = save_pars(all = TRUE)
  )
  print(fit3)
  me3 <- conditional_effects(fit3, method = "predict")
  expect_ggplot(plot(me3, ask = FALSE)[[2]])
  expect_range(LOO(fit3)$estimates[3, 1], 650, 740)

  # posterior checks of the censored model
  expect_ggplot(SW(pp_check(
    fit3, type = 'dens_overlay_grouped', group = 'sex', ndraws = 10
  )))
  expect_ggplot(SW(pp_check(
    fit3, type = 'intervals', x = 'patient', ndraws = NULL
  )))
  expect_ggplot(SW(pp_check(fit3, type = 'loo_intervals', ndraws = NULL)))
  expect_ggplot(SW(pp_check(fit3, type = 'loo_pit_overlay', ndraws = 10)))

  # enables rstan specific functionality
  fit3 <- add_rstan_model(fit3)
  expect_range(LOO(fit3, moment_match = TRUE)$estimates[3, 1], 650, 740)
  bridge <- bridge_sampler(fit3)
  expect_true(is.numeric(bridge$logml))

  testthat::skip_if_not_installed("ggfortify")
  # TODO: automatically extract status_y within pp_check
  expect_ggplot(SW(pp_check(
    fit3, type = 'km_overlay', ndraws = 10,
    status_y = 1 - kidney$censored
  )))
}))

test_that("Binomial model from brm doc works correctly", suppressWarnings({
  ## Probit regression using the binomial family
  n <- sample(1:10, 100, TRUE)  # number of trials
  success <- rbinom(100, size = n, prob = 0.4)
  x <- rnorm(100)
  data4 <- data.frame(n, success, x)
  fit4 <- brm(
    success | trials(n) ~ x, data = data4,
    family = binomial("probit"), refresh = 0
  )
  print(fit4)
  ce <- conditional_effects(fit4)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
}))

test_that("Non-linear model from brm doc works correctly", suppressWarnings({
  x <- abs(rnorm(100))
  y <- rnorm(100, mean = 2 - 1.5^x, sd = 1)
  data5 <- data.frame(x, y)
  fit5 <- brm(
    bf(y ~ a1 - a2^x, a1 + a2 ~ 1, nl = TRUE), data = data5,
    prior = prior(normal(0, 2), nlpar = a1) +
      prior(normal(0, 2), nlpar = a2), refresh = 0
  )
  print(fit5)
  ce <- conditional_effects(fit5)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
}))

test_that("Models from hypothesis doc work correctly", suppressWarnings({
  prior <- c(set_prior("normal(0,2)", class = "b"),
             set_prior("student_t(10,0,1)", class = "sigma"),
             set_prior("student_t(10,0,1)", class = "sd"))

  ## fit a linear mixed effects models
  fit <- brm(time ~ age + sex + disease + (1 + age|patient),
             data = kidney, family = lognormal(),
             prior = prior, sample_prior = TRUE,
             control = list(adapt_delta = 0.95),
             refresh = 0)
  print(fit)

  ## perform two-sided hypothesis testing
  hyp1 <- hypothesis(fit, "sexfemale = age + diseasePKD")
  expect_range(hyp1$hypothesis$Estimate, 0.6, 0.74)
  expect_range(hyp1$hypothesis$Evid.Ratio, 1, 5)
  expect_true(is(plot(hyp1)[[1]], "ggplot"))

  ## perform one-sided hypothesis testing
  hyp2 <- hypothesis(fit, "diseasePKD + diseaseGN - 3 < 0")
  expect_range(hyp2$hypothesis$Evid.Ratio, 300)

  ## test more than one hypothesis at once
  hyp3 <- c("diseaseGN = diseaseAN", "2 * diseaseGN - diseasePKD = 0")
  hyp3 <- hypothesis(fit, hyp3)
  expect_equal(dim(hyp3$hypothesis), c(2, 8))
}))

test_that("bridgesampling methods work correctly", suppressWarnings({
  # model with the treatment effect
  fit1 <- brm(
    count ~ zAge + zBase + Trt,
    data = epilepsy, family = negbinomial(),
    prior = prior(normal(0, 1), class = b),
    save_pars = save_pars(all = TRUE), refresh = 0,
    backend = "rstan"
  )
  print(fit1)
  # model without the treatment effect
  fit2 <- brm(
    count ~ zAge + zBase,
    data = epilepsy, family = negbinomial(),
    prior = prior(normal(0, 1), class = b),
    save_pars = save_pars(all = TRUE), refresh = 0,
    backend = "rstan"
  )
  print(fit2)

  # compute the bayes factor
  expect_gt(bayes_factor(fit2, fit1)$bf, 1)
  # compute the posterior model probabilities
  pp1 <- post_prob(fit2, fit1)
  expect_gt(pp1[1], pp1[2])
  # specify prior model probabilities
  pp2 <- post_prob(fit2, fit1, prior_prob = c(0.8, 0.2))
  expect_gt(pp2[1], pp1[1])
}))

test_that("varying slopes without overall effects work", suppressWarnings({
  epilepsy$visit_num <- as.numeric(epilepsy$visit)
  fit1 <- brm(count ~ zAge + zBase * Trt + (visit_num | patient),
              data = epilepsy, family = gaussian(),
              chains = 2, refresh = 0)
  print(fit1)
  ce <- conditional_effects(fit1)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  # test reloo
  loo1 <- LOO(fit1)
  reloo1 <- reloo(loo1, fit1, chains = 1, iter = 100)
  expect_range(reloo1$estimates[3, 1], 1600, 1700)
  reloo2 <- LOO(fit1, reloo = TRUE, chains = 1, iter = 100)
  expect_range(reloo2$estimates[3, 1], 1600, 1700)

  conditions <- data.frame(zAge = 0, zBase = 0, Trt = 0)
  ce <- conditional_effects(fit1, conditions = conditions)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])

  conditions <- data.frame(zAge = 0, zBase = 0, patient = 0,
                           Trt = 0, visit_num = 1:5)
  ce <- conditional_effects(fit1, conditions = conditions, re_formula = NULL)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])

  expect_range(WAIC(fit1)$estimates[3, 1], 1500, 1600)
  expect_equal(dim(predict(fit1)), c(nobs(fit1), 4))
}))

test_that("update works correctly for some special cases", suppressWarnings({
  # models are recompiled when changing number of FEs from 0 to 1
  fit1 <- brm(count ~ 1, data = epilepsy, refresh = 0)
  fit2 <- update(fit1, ~ . + Trt, newdata = epilepsy, refresh = 0)
  expect_equal(rownames(fixef(fit2)), c("Intercept", "Trt1"))
  fit3 <- update(fit2, ~ . - Trt, newdata = epilepsy, refresh = 0)
  expect_equal(rownames(fixef(fit3)), c("Intercept"))

  # test that family is correctly updated
  fit4 <- update(fit1, family = student(), refresh = 0)
  expect_equal(family(fit4)$family, "student")
  expect_equal(formula(fit4)$family$family, "student")
  fit5 <- update(fit1, bf(~., family = student()), refresh = 0)
  expect_equal(family(fit5)$family, "student")
  expect_equal(formula(fit5)$family$family, "student")
}))

test_that("Fixing parameters to constants works correctly", suppressWarnings({
  bprior <- prior(normal(0, 1), class = "b") +
    prior(constant(2), class = "b", coef = "zBase") +
    prior(constant(0.5), class = "sd")

  fit <- brm(count ~ zAge + zBase + (1 | patient),
             data = epilepsy, prior = bprior)
  print(summary(fit))
  expect_range(waic(fit)$estimates[3, 1], 1790, 1840)
  ce <- conditional_effects(fit)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
}))
