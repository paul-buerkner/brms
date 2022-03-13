source("setup_local_tests.R")

test_that("Poisson model from brm doc works correctly", {
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
})

test_that("Ordinal model from brm doc works correctly", {
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
  expect_warning(ce <- conditional_effects(fit2, effect = "treat"),
                 "Predictions are treated as continuous variables")
  expect_ggplot(plot(ce)[[1]])
})

test_that("Survival model from brm doc works correctly", {
  ## Survival regression modeling the time between the first
  ## and second recurrence of an infection in kidney patients.
  fit3 <- brm(
    time | cens(censored) ~ age * sex + disease + (1|patient),
    data = kidney, family = lognormal(), refresh = 0,
    threads = threading(2, grainsize = 100),
    backend = "cmdstanr"
  )
  print(fit3)
  me3 <- conditional_effects(fit3, method = "predict")
  expect_ggplot(plot(me3, ask = FALSE)[[2]])
  expect_range(LOO(fit3)$estimates[3, 1], 650, 740)

  # enables rstan specific functionality
  fit3 <- add_rstan_model(fit3)
  expect_range(LOO(fit3, moment_match = TRUE)$estimates[3, 1], 650, 740)
  bridge <- bridge_sampler(fit3)
  expect_true(is.numeric(bridge$logml))
})

test_that("Binomial model from brm doc works correctly", {
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
})

test_that("Non-linear model from brm doc works correctly", {
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
})

test_that("ARMA models work correctly", {
  set.seed(1234)
  N <- 100
  y <- arima.sim(list(ar = c(0.7, -0.5, 0.04, 0.2, -0.4)), N)
  dat <- list(y = y, x = rnorm(N), g = sample(1:5, N, TRUE))

  fit_ar <- brm(
    y ~ x + ar(p = 5), data = dat,
    prior = prior(normal(0, 5), class = "ar"),
    chains = 2, refresh = 0
  )
  print(fit_ar)
  ar <- colMeans(as.matrix(fit_ar, variable = "ar"))
  expect_range(ar[1], 0.5, 0.9)
  expect_range(ar[3], -0.2, 0.25)
  expect_range(ar[5], -0.6, -0.1)
  expect_ggplot(plot(conditional_effects(fit_ar))[[1]])

  fit_ma <- brm(y ~ x + ma(q = 1), data = dat,
                chains = 2, refresh = 0)
  print(fit_ma)
  expect_gt(LOO(fit_ma)$estimates[3, 1], LOO(fit_ar)$estimates[3, 1])

  fit_arma <- brm(
    y ~ x + (1|g) + arma(gr = g, cov = TRUE), data = dat,
    prior = prior(normal(0, 5), class = "ar") +
      prior(normal(0, 6), class = "ma"),
    chains = 2, refresh = 0
  )
  print(fit_arma)
  expect_range(waic(fit_arma)$estimates[3, 1], 250, 350)
  expect_equal(dim(predict(fit_arma)), c(nobs(fit_arma), 4))
  expect_ggplot(plot(conditional_effects(fit_arma), plot = FALSE)[[1]])

  fit_arma_pois <- brm(
    count ~ Trt + (1 | patient) + arma(visit, patient, cov = TRUE),
    data = epilepsy, family = poisson(),
    chains = 2, refresh = 0
  )
  print(fit_arma_pois)
  expect_range(waic(fit_arma_pois)$estimates[3, 1], 1100, 1200)
  expect_equal(dim(predict(fit_arma_pois)), c(nobs(fit_arma_pois), 4))
  expect_ggplot(plot(conditional_effects(fit_arma_pois), plot = FALSE)[[1]])
})

test_that("Models from hypothesis doc work correctly", {
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
})

test_that("categorical models work correctly", {
  fit2 <- brm(rating ~ period + carry + treat + (1|test|subject),
              data = inhaler, family = categorical, iter = 500,
              prior = c(prior(normal(0,5), "b"),
                        prior(normal(0,5), "Intercept")),
              chains = 2, refresh = 0)
  print(fit2)
  expect_range(WAIC(fit2)$estimates[3, 1], 830, 900)
  ncat <- length(unique(inhaler$rating))
  expect_equal(dim(predict(fit2)), c(nobs(fit2), ncat))
  expect_equal(dim(fitted(fit2)), c(nobs(fit2), 4, ncat))
  expect_equal(dim(fitted(fit2, scale = "linear")),
               c(nobs(fit2), 4, ncat - 1))

  # tests with new data
  newd <- inhaler[1:10, ]
  newd$rating <- NULL
  expect_equal(dim(predict(fit2, newdata = newd)), c(10, ncat))

  ce <- conditional_effects(fit2, categorical = TRUE)
  expect_ggplot(plot(ce, plot = FALSE)[[1]])
})

test_that("bridgesampling methods work correctly", {
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
})

test_that("varying slopes without overall effects work", {
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

  conditions <- data.frame(zAge = 0, zBase = 0,
                           Trt = 0, visit = c(1:4, NA))
  ce <- conditional_effects(fit1, conditions = conditions, re_formula = NULL)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])

  expect_range(WAIC(fit1)$estimates[3, 1], 1500, 1600)
  expect_equal(dim(predict(fit1)), c(nobs(fit1), 4))
})

test_that("multivariate normal models work correctly", {
  set.seed(1234)
  N <- 300
  y1 <- rnorm(N)
  y2 <- rnorm(N, mean = 1, sd = 2)
  x <- rnorm(N, 1)
  month <- sample(1:36, N, replace = TRUE)
  id <- sample(1:10, N, replace = TRUE)
  tim <- sample(1:N, N)
  data <- data.frame(y1, y2, x, month, id, tim)

  fit_mv1 <- brm(mvbind(y1, y2) ~ s(x) + poly(month, 3) +
                   (1|x|id) + arma(tim, id, q = 0),
                 data = data,
                 prior = c(prior_(~normal(0,5), resp = "y1"),
                           prior_(~normal(0,5), resp = "y2"),
                           prior_(~lkj(5), class = "rescor")),
                 sample_prior = TRUE,
                 iter = 1000, chains = 2, refresh = 0)
  print(fit_mv1)

  expect_equal(dim(predict(fit_mv1)), c(300, 4, 2))
  expect_equal(dim(fitted(fit_mv1)), c(300, 4, 2))
  newdata <- data.frame(month = 1, y1 = 0, y2 = 0, x = 0, id = 1, tim = 1)
  expect_equal(dim(predict(fit_mv1, newdata = newdata)), c(1, 4, 2))
  cs <- conditional_smooths(fit_mv1, ndraws = 750)
  expect_equal(length(cs), 2)
  expect_ggplot(plot(cs, ask = FALSE)[[2]])

  fit_mv2 <- brm(mvbind(y1, y2) ~ 1, data = data,
                 prior = prior_(~lkj(5), class = "rescor"),
                 sample_prior = TRUE, iter = 1000, refresh = 0)
  print(fit_mv2)
  waic_mv <- WAIC(fit_mv1, fit_mv2, ndraws = 100)
  expect_true(waic_mv$ic_diffs__[1, "WAIC"] > 0)
})

test_that("generalized multivariate models work correctly", {
  data("BTdata", package = "MCMCglmm")
  bform <- (bf(tarsus ~ sex + (1|p|fosternest)) + skew_normal()) +
    (bf(back ~ s(tarsus, by = sex) + (1|p|fosternest)) + gaussian())
  fit_mv <- brm(bform, BTdata, chains = 2, iter = 1000, refresh = 0)

  print(fit_mv)
  expect_ggplot(pp_check(fit_mv, resp = "back"))
  expect_range(waic(fit_mv)$estimates[3, 1], 4300, 4400)
  expect_ggplot(plot(conditional_effects(fit_mv), ask = FALSE)[[1]])
  expect_ggplot(plot(conditional_smooths(fit_mv))[[1]])
  expect_equal(dim(coef(fit_mv)$fosternest), c(104, 4, 7))
})

test_that("ZI and HU models work correctly", {
  fit_hu <- brm(
    bf(count ~ zAge + zBase * Trt + (1|id1|patient),
       hu ~ zAge + zBase * Trt + (1|id1|patient)),
    data = epilepsy, family = hurdle_poisson(),
    prior = c(prior(normal(0, 5)),
              prior(normal(0, 5), dpar = "hu")),
    chains = 2, refresh = 0
  )
  print(fit_hu)
  expect_equal(dim(predict(fit_hu)), c(nobs(fit_hu), 4))
  expect_ggplot(plot(conditional_effects(fit_hu), ask = FALSE)[[2]])

  fit_zi <- brm(
    bf(count ~ zAge + zBase * Trt + (1|patient), zi ~ Trt),
    data = epilepsy, family = zero_inflated_negbinomial(),
    prior = prior(normal(0,5)) +
      prior(normal(0,3), class = "sd") +
      prior(cauchy(0,5), class = "shape") +
      prior(normal(0,5), dpar = "zi") ,
    chains = 2, refresh = 0
  )
  print(fit_zi)
  expect_equal(dim(predict(fit_zi)), c(nobs(fit_zi), 4))
  expect_ggplot(plot(conditional_effects(fit_zi), ask = FALSE)[[2]])
  waic_zi <- WAIC(fit_hu, fit_zi, ndraws = 100)
  expect_equal(dim(waic_zi$ic_diffs__), c(1, 2))

  ## zero_inflated beta model
  data("GasolineYield", package = "betareg")
  dat <- GasolineYield
  dat$yield[c(1, 5, 8, 12, 16)] <- 0
  fit_zibeta <- brm(
    yield ~ batch + temp, data = dat,
    family = zero_inflated_beta(),
    chains = 2, inits = 0, refresh = 0
  )
  print(fit_zibeta)
  expect_equal(dim(predict(fit_zibeta)), c(nobs(fit_zibeta), 4))
  expect_ggplot(plot(conditional_effects(fit_zibeta), ask = FALSE)[[1]])
  expect_range(WAIC(fit_zibeta)$estimates[3, 1], -100, -70)
})

test_that("Non-linear models work correctly", {
  fit_loss <- brm(
    bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
       ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
       nl = TRUE),
    data = loss, family = gaussian(),
    prior = c(prior(normal(5000, 1000), nlpar = "ult"),
              prior(normal(1, 2), nlpar = "omega"),
              prior(normal(45, 10), nlpar = "theta")),
    control = list(adapt_delta = 0.9),
    refresh = 0
  )
  print(fit_loss)
  expect_ggplot(plot(conditional_effects(fit_loss))[[1]])
  expect_range(LOO(fit_loss)$estimates[3, 1], 700, 720)
})

test_that("Non-linear models of distributional parameters work correctly", {
  set.seed(1234)
  x <- rnorm(100)
  y <- rnorm(100, 1, exp(x))
  dat <- data.frame(y, x, g = rep(1:10, each = 10))
  bform <- bf(y ~ x + (1|V|g)) +
    nlf(sigma ~ a) +
    lf(a ~ x + (1|V|g)) +
    gaussian()
  bprior <- prior(normal(0, 3), nlpar = a)
  fit <- brm(bform, dat, prior = bprior, chains = 2, refresh = 0)
  print(fit)
  ce <- conditional_effects(fit, method = "predict")
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_equal(dim(fitted(fit, dat[1:10, ])), c(10, 4))
  expect_range(LOO(fit)$estimates[3, 1], 240, 350)
})

test_that("Nested non-linear models work correctly", {
  set.seed(2345)
  dat <- data.frame(x = rnorm(300))
  dat$y <- 0.3 + 0.7 * brms:::inv_logit(2 * dat$x)
  bform <- bf(
    y ~ lb + (1 - lb) * inv_logit(b * x),
    b + a ~ 1, nlf(lb ~ inv_logit(a)),
    nl = TRUE
  )
  bprior <- prior(normal(0, 1), nlpar = "a") +
    prior(normal(0, 1), nlpar = "b")
  fit <- brm(bform, dat, prior = bprior, family = Beta(), refresh = 0)
  print(fit)

  ce <- conditional_effects(fit)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(bayes_R2(fit)[, 1], 0.2, 0.55)
})

test_that("Multivariate GAMMs work correctly", {
  set.seed(4312)
  n <- 200
  sig <- 2
  dat <- mgcv::gamSim(1, n = n, scale = sig)
  fit_gam <- brm(
    y ~ t2(x0, x2) + s(x1), data = dat, chains = 2,
    control = list(adapt_delta = 0.95), refresh = 0
  )
  print(fit_gam)

  ce <- conditional_effects(fit_gam)
  expect_ggplot(plot(ce, ask = FALSE, rug = TRUE, points = TRUE)[[1]])
  ce <- conditional_effects(fit_gam, surface = TRUE, too_far = 0.05)
  expect_ggplot(plot(ce, ask = FALSE, rug = TRUE)[[1]])

  cs <- conditional_smooths(fit_gam, resolution = 25)
  expect_ggplot(plot(cs, rug = TRUE, ask = FALSE)[[1]])
  cs <- conditional_smooths(fit_gam, resolution = 100, too_far = 0.05)
  expect_ggplot(plot(cs, rug = TRUE, ask = FALSE)[[1]])

  expect_range(loo(fit_gam)$estimates[3, 1], 830, 870)
  expect_equal(dim(predict(fit_gam)), c(nobs(fit_gam), 4))

  newd <- data.frame(x0 = (0:30)/30, x1 = (0:30)/30,
                     x2 = (0:30)/30, x3 = (0:30)/30)
  prfi <- cbind(predict(fit_gam, newd), fitted(fit_gam, newdata = newd))
  expect_range(prfi[, 1], prfi[, 5] - 0.25, prfi[, 5] + 0.25)
})

test_that("GAMMs with factor variable in 'by' work correctly", {
  set.seed(7)
  dat <- mgcv::gamSim(4, n = 200, dist = "normal")
  fit_gam2 <- brm(y ~ fac + s(x2, by = fac, k = 4), dat,
                  chains = 2, refresh = 0)
  print(fit_gam2)

  ce <- conditional_effects(fit_gam2, "x2:fac")
  expect_ggplot(plot(ce, points = TRUE, ask = FALSE)[[1]])
  cs <- conditional_smooths(fit_gam2, res = 10)
  expect_ggplot(plot(cs, rug = TRUE, ask = FALSE)[[1]])

  fit_gam3 <- brm(y ~ fac + t2(x1, x2, by = fac), dat,
                  chains = 2, refresh = 0)
  print(fit_gam3)

  ce <- conditional_effects(fit_gam3, "x2:fac")
  expect_ggplot(plot(ce, points = TRUE, ask = FALSE)[[1]])
  cs <- conditional_smooths(fit_gam3, too_far = 0.1)
  expect_ggplot(plot(cs, rug = TRUE, ask = FALSE)[[1]])
  expect_ggplot(plot(cs, rug = TRUE, stype = "raster", ask = FALSE)[[1]])
})

test_that("generalized extreme value models work correctly", {
  data(fremantle, package = "ismev")
  fremantle <- transform(fremantle, cYear = Year - median(Year))
  knots <- with(
    fremantle,
    list(cYear = c(
      min(Year) - c(10, 0), 1945,
      max(Year) + c(0, 10)) - median(Year)
    )
  )

  fit_gev <- brm(
    bf(SeaLevel ~ cYear + SOI,
       sigma ~ s(cYear, bs = "bs", m = 1, k = 3) + SOI),
    data = fremantle, family = gen_extreme_value(),
    knots = knots, inits = 0.5, chains = 4,
    control = list(adapt_delta = 0.95), refresh = 0
  )
  print(fit_gev)

  prfi <- cbind(predict(fit_gev), fitted(fit_gev))
  expect_range(prfi[, 1], prfi[, 5] - 0.03, prfi[, 5] + 0.03)
  # expect_range(loo(fit_gev)$estimates[3, 1], -115, -95)
  ce <- conditional_effects(fit_gev, "cYear")
  expect_ggplot(plot(ce, points = TRUE, ask = FALSE)[[1]])
})

test_that("update works correctly for some special cases", {
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
})

test_that("Wiener diffusion models work correctly", {
  set.seed(312)
  x <- rnorm(100, mean = 1)
  dat <- rwiener(n = 1, alpha = 2, tau = .3, beta = .5, delta = .5 + x)
  dat$x <- x

  fit_d1 <- brm(bf(q | dec(resp) ~ x), dat,
                family = wiener(), refresh = 0)
  print(fit_d1)
  expect_ggplot(plot(conditional_effects(fit_d1), ask = FALSE)[[1]])
  expect_ggplot(pp_check(fit_d1))
  pp <- posterior_predict(fit_d1, ndraws = 10, negative_rt = TRUE)
  expect_true(min(pp) < 0)

  fit_d2 <- brm(bf(q | dec(resp) ~ x, ndt ~ x),
                dat, family = wiener(), refresh = 0)
  print(fit_d2)
  expect_ggplot(plot(conditional_effects(fit_d2), ask = FALSE)[[1]])
  expect_ggplot(pp_check(fit_d2))
  pred <- predict(fit_d2, ndraws = 10, negative_rt = TRUE, summary = FALSE)
  expect_true(any(pred < 0))

  waic_d <- WAIC(fit_d1, fit_d2)
  expect_equal(dim(waic_d$ic_diffs__), c(1, 2))
})

test_that("disc parameter in ordinal models is handled correctly", {
  fit <- brm(
    bf(rating ~ period + carry + treat + (1|subject), disc ~ 1),
    data = inhaler, family = cumulative(),
    prior = prior(normal(0,5)),
    chains = 2, refresh = 0
  )
  print(fit)
  expect_range(waic(fit)$estimates[3, 1], 870, 920)
  ncat <- length(unique(inhaler$rating))
  expect_equal(dim(predict(fit)), c(nobs(fit), ncat))
  expect_ggplot(plot(
    conditional_effects(fit), ask = FALSE,
    points = TRUE, point_args = list(width = 0.3)
  )[[3]])
})

test_that("Argument `incl_thres` works correctly for non-grouped thresholds", {
  fit <- brm(
    bf(rating ~ period + carry + treat + (1|subject)),
    data = inhaler, family = cumulative(),
    prior = prior(normal(0,5)),
    chains = 2, refresh = 0
  )
  thres_minus_eta <- posterior_linpred(fit, incl_thres = TRUE)
  bprep <- prepare_predictions(fit)
  thres <- bprep$thres$thres
  eta <- posterior_linpred(fit)
  thres_minus_eta_ch <- apply(thres, 2, "-", eta)
  thres_minus_eta_ch <- array(
    thres_minus_eta_ch, dim = c(nrow(thres), ncol(eta), ncol(thres))
  )
  expect_equivalent(thres_minus_eta, thres_minus_eta_ch)
})

test_that("Mixture models work correctly", {
  set.seed(12346)
  dat <- data.frame(
    y = c(rnorm(300), rnorm(100, 6), rnorm(200, 12)),
    x = rnorm(600),
    z = sample(0:1, 600, TRUE)
  )
  bform1 <- bf(
    y ~ 1, mu1 + mu2 ~ x, mu3 ~ z,
    sigma2 = "sigma1", sigma3 = "sigma1"
  )
  mixfam <- mixture(gaussian(), nmix = 3)
  prior <- c(
    prior(normal(0, 5), Intercept, dpar = mu1),
    prior(normal(5, 5), Intercept, dpar = mu2),
    prior(normal(10, 5), Intercept, dpar = mu3)
  )

  fit1 <- brm(
    bform1, data = dat, family = mixfam,
    prior = c(prior, prior(dirichlet(1, 1, 1), theta)),
    chains = 2, inits = 0, refresh = 0
  )
  print(fit1)
  expect_ggplot(pp_check(fit1))
  loo1 <- LOO(fit1)
  expect_equal(dim(pp_mixture(fit1)), c(nobs(fit1), 4, 3))

  bform2 <- bf(bform1, theta1 = 1, theta2 = 1, theta3 = 1)
  fit2 <- brm(
    bform2, data = dat, family = mixfam,
    prior = prior, chains = 2,
    inits = 0, refresh = 0
  )
  print(fit2)
  expect_ggplot(pp_check(fit2))
  loo2 <- LOO(fit2)
  expect_gt(loo2$estimates[3, 1], loo1$estimates[3, 1])
  expect_equal(dim(pp_mixture(fit2)), c(nobs(fit2), 4, 3))

  # MCMC chains get stuck when fitting this model
  # bform3 <- bf(bform1, theta1 ~ z, theta2 ~ 1)
  # prior3 <- prior +
  #   prior(normal(0, 1), dpar = theta1) +
  #   prior(normal(0, 1), Intercept, dpar = theta1) +
  #   prior(normal(0, 1), Intercept, dpar = theta2)
  # fit3 <- brm(
  #   bform3, data = dat, family = mixfam,
  #   prior = prior3, init_r = 0.1,
  #   chains = 1, refresh = 0
  # )
  # print(fit3)
  # expect_ggplot(pp_check(fit3))
  # loo3 <- LOO(fit3, pointwise = TRUE)
  # expect_range(loo3$estimates[3, 1],
  #   loo1$estimates[3, 1] - 20, loo1$estimates[3, 1] + 20
  # )
  # expect_equal(dim(pp_mixture(fit3)), c(nobs(fit3), 4, 3))
})

test_that("Gaussian processes work correctly", {
  ## Basic GPs
  set.seed(1112)
  dat <- mgcv::gamSim(1, n = 30, scale = 2)
  fit1 <- brm(y ~ gp(x0) + x1 + gp(x2) + x3, dat,
              chains = 2, refresh = 0)
  print(fit1)
  expect_ggplot(pp_check(fit1))
  ce <- conditional_effects(fit1, ndraws = 200, nug = 1e-07)
  expect_ggplot(plot(ce, ask = FALSE)[[3]])
  expect_range(WAIC(fit1)$estimates[3, 1], 100, 200)

  # multivariate GPs
  fit2 <- brm(y ~ gp(x1, x2), dat, chains = 2, refresh = 0)
  print(fit2)
  expect_ggplot(pp_check(fit2))
  ce <- conditional_effects(
    fit2, ndraws = 200, nug = 1e-07,
    surface = TRUE, resolution = 10
  )
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit2)$estimates[3, 1], 100, 200)

  # GP with continuous 'by' variable
  fit3 <- brm(y ~ gp(x1, by = x2), dat, chains = 2, refresh = 0)
  print(fit3)
  expect_ggplot(pp_check(fit3))
  ce <- conditional_effects(fit3, ndraws = 200, nug = 1e-07)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit3)$estimates[3, 1], 100, 200)

  # GP with factor 'by' variable
  dat2 <- mgcv::gamSim(4, n = 100, scale = 2)
  fit4 <- brm(y ~ gp(x2, by = fac), dat2, chains = 2, refresh = 0)
  print(fit4)
  expect_ggplot(pp_check(fit4))
  ce <- conditional_effects(fit4, ndraws = 200, nug = 1e-07)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit4)$estimates[3, 1], 400, 600)
})

test_that("Approximate Gaussian processes work correctly", {
  set.seed(1245)
  dat <- mgcv::gamSim(4, n = 200, scale = 2)

  # isotropic approximate GP
  fit1 <- brm(
    y ~ gp(x1, x2, by = fac, k = 10, c = 5/4),
    data = dat, chains = 2, cores = 2, refresh = 0,
    control = list(adapt_delta = 0.99)
  )
  print(fit1)
  expect_range(bayes_R2(fit1)[1, 1], 0.45, 0.7)
  ce <- conditional_effects(
    fit1, "x2:x1", conditions = data.frame(fac = unique(dat$fac)),
    resolution = 20, surface = TRUE
  )
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit1)$estimates[3, 1], 900, 1000)

  # non isotropic approximate GP
  fit2 <- brm(
    y ~ gp(x1, x2, by = fac, k = 10, c = 5/4, iso = FALSE),
    data = dat, chains = 2, cores = 2, refresh = 0,
    control = list(adapt_delta = 0.99)
  )
  print(fit2)
  expect_range(bayes_R2(fit2)[1, 1], 0.50, 0.62)
  ce <- conditional_effects(
    fit2, "x2:x1", conditions = data.frame(fac = unique(dat$fac)),
    resolution = 20, surface = TRUE
  )
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit2)$estimates[3, 1], 870, 970)
})

test_that("SAR models work correctly", {
  data(oldcol, package = "spdep")
  fit_lagsar <- brm(CRIME ~ INC + HOVAL + sar(COL.nb),
                    data = COL.OLD, data2 = list(COL.nb = COL.nb),
                    chains = 2, refresh = 0)
  print(fit_lagsar)
  expect_ggplot(pp_check(fit_lagsar))
  ce <- conditional_effects(fit_lagsar, ndraws = 200)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit_lagsar)$estimates[3, 1], 350, 380)

  fit_errorsar <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "error"),
                      data = COL.OLD, data2 = list(COL.nb = COL.nb),
                      chains = 2, refresh = 0)
  print(fit_errorsar)
  expect_ggplot(pp_check(fit_errorsar))
  ce <- conditional_effects(fit_errorsar, ndraws = 200)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit_errorsar)$estimates[3, 1], 350, 380)
})

test_that("CAR models work correctly", {
  # generate some spatial data
  set.seed(4331)
  east <- north <- 1:10
  Grid <- expand.grid(east, north)
  K <- nrow(Grid)

  # set up distance and neighbourhood matrices
  distance <- as.matrix(dist(Grid))
  W <- array(0, c(K, K))
  W[distance == 1] <- 1
  rownames(W) <- 1:nrow(W)

  # generate the covariates and response data
  x1 <- rnorm(K)
  x2 <- rnorm(K)
  theta <- rnorm(K, sd = 0.05)
  phi <- rmulti_normal(
    1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
  )
  eta <- x1 + x2 + phi
  prob <- exp(eta) / (1 + exp(eta))
  size <- rep(50, K)
  y <- rbinom(n = K, size = size, prob = prob)
  dat <- data.frame(y, size, x1, x2, obs = 1:length(y))

  # fit a CAR model
  fit_car <- brm(
    y | trials(size) ~ x1 + x2 + car(W, obs),
    data = dat, data2 = list(W = W), family = binomial(),
    chains = 2, refresh = 0
  )
  print(fit_car)
  expect_ggplot(pp_check(fit_car))
  ce = conditional_effects(fit_car, ndraws = 200)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(LOO(fit_car)$estimates[3, 1], 450, 550)
  expect_false(isTRUE(all.equal(
    fitted(fit_car, newdata = dat[1:5, ]),
    fitted(fit_car, newdata = dat[1:5, ], incl_autocor = FALSE)
  )))

  newdata <- data.frame(x1 = 0, x2 = 0, size = 50, obs = 1)
  pp <- posterior_predict(fit_car, newdata = newdata)
  expect_equal(dim(pp), c(ndraws(fit_car), 1))

  newdata <- data.frame(x1 = 0, x2 = 0, size = 50, obs = -1)
  new_W <- W
  rownames(W)[1] <- "-1"
  newdata2 <- list(W = new_W)
  expect_error(predict(fit_car, newdata = newdata, newdata2 = newdata2),
               "Cannot handle new locations in CAR models")
})

test_that("Missing value imputation works correctly", {
  library(mice)
  data("nhanes", package = "mice")

  # missing value imputation via multiple imputation
  imp <- mice(nhanes)
  fit_imp1 <- brm_multiple(bmi ~ age * chl, imp, chains = 1,
                           backend = "rstan", refresh = 0)
  print(fit_imp1)
  expect_equal(ndraws(fit_imp1), 5000)
  expect_equal(dim(fit_imp1$rhats), c(5, length(variables(fit_imp1))))

  fit_imp1 <- update(fit_imp1, . ~ chl, newdata = imp)
  print(fit_imp1)
  expect_true(!"b_age" %in% variables(fit_imp1))
  expect_equal(ndraws(fit_imp1), 5000)

  # missing value imputation within Stan
  bform <- bf(bmi | mi() ~ age * mi(chl)) +
    bf(chl | mi() ~ age) + set_rescor(FALSE)
  fit_imp2 <- brm(bform, data = nhanes, backend = "rstan", refresh = 0)
  print(fit_imp2)
  pred <- predict(fit_imp2)
  expect_true(!anyNA(pred))
  ce <- conditional_effects(fit_imp2, resp = "bmi")
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  loo <- LOO(fit_imp2, newdata = na.omit(fit_imp2$data))
  expect_range(loo$estimates[3, 1], 200, 220)

  # overimputation within Stan
  dat <- nhanes
  dat$sdy <- 5
  bform <- bf(bmi | mi() ~ age * mi(chl)) +
    bf(chl | mi(sdy) ~ age) + set_rescor(FALSE)
  fit_imp3 <- brm(bform, data = dat,
                  save_pars = save_pars(latent = TRUE),
                  backend = "rstan", refresh = 0)
  print(fit_imp3)
  pred <- predict(fit_imp3)
  expect_true(!anyNA(pred))
  ce <- conditional_effects(fit_imp3, resp = "bmi")
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  loo <- LOO(fit_imp3, newdata = na.omit(fit_imp3$data))
  expect_range(loo$estimates[3, 1], 200, 220)
})

test_that("student-t-distributed group-level effects work correctly", {
  fit <- brm(
    count ~ Trt * zBase + (1 | gr(patient, dist = "student")),
    data = epilepsy, family = poisson(),
    chains = 1, refresh = 0
  )
  print(summary(fit))
  expect_true("df_patient" %in% variables(fit))
  expect_true(!"udf_1" %in% variables(fit))
  waic <- suppressWarnings(waic(fit))
  expect_range(waic$estimates[3, 1], 1300, 1400)
})

test_that("multinomial models work correctly", {
  set.seed(1245)
  N <- 100
  dat <- data.frame(
    y1 = rbinom(N, 10, 0.1), y2 = rbinom(N, 10, 0.4),
    y3 = rbinom(N, 10, 0.7), x = rnorm(N)
  )
  dat$size <- with(dat, y1 + y2 + y3)
  dat$y <- with(dat, cbind(y1, y2, y3))

  fit <- brm(y | trials(size) ~ x, data = dat,
             family = multinomial(), refresh = 0)
  print(summary(fit))
  pred <- predict(fit)
  expect_equal(dim(pred), c(nobs(fit), 4, 3))
  expect_equal(dimnames(pred)[[3]], c("y1", "y2", "y3"))
  waic <- waic(fit)
  expect_range(waic$estimates[3, 1], 550, 600)
  ce <- conditional_effects(fit, categorical = TRUE)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
})

test_that("dirichlet models work correctly", {
  set.seed(1246)
  N <- 100
  dat <- as.data.frame(rdirichlet(N, c(10, 5, 1)))
  names(dat) <- c("y1", "y2", "y3")
  dat$x <- rnorm(N)
  dat$y <- with(dat, cbind(y1, y2, y3))

  fit <- brm(y ~ x, data = dat, family = dirichlet(), refresh = 0)
  print(summary(fit))
  expect_output(print(fit), "muy2 = logit")
  pred <- predict(fit)
  expect_equal(dim(pred), c(nobs(fit), 4, 3))
  expect_equal(dimnames(pred)[[3]], c("y1", "y2", "y3"))
  waic <- waic(fit)
  expect_range(waic$estimates[3, 1], -530, -500)
  ce <- conditional_effects(fit, categorical = TRUE)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
})

test_that("logistic_normal models work correctly", {
  set.seed(1246)
  N <- 100
  dat <- as.data.frame(rdirichlet(N, c(10, 5, 1)))
  names(dat) <- c("y1", "y2", "y3")
  dat$x <- rnorm(N)
  dat$y <- with(dat, cbind(y1, y2, y3))

  fit <- brm(y ~ x, data = dat, family = logistic_normal(), refresh = 0)
  print(summary(fit))
  expect_output(print(fit), "muy2 = identity")
  pred <- predict(fit, ndraws = 250)
  expect_equal(dim(pred), c(nobs(fit), 4, 3))
  expect_equal(dimnames(pred)[[3]], c("y1", "y2", "y3"))
  waic <- waic(fit, ndraws = 250)
  expect_range(waic$estimates[3, 1], -530, -460)
})

test_that("Addition argument 'subset' works correctly", {
  set.seed(12454)
  data("BTdata", package = "MCMCglmm")
  BTdata$sub1 <- sample(0:1, nrow(BTdata), replace = TRUE)
  BTdata$sub2 <- sample(0:1, nrow(BTdata), replace = TRUE)

  bform <- bf(tarsus | subset(sub1) ~ sex + (1|p|fosternest) + (1|q|dam)) +
    bf(back | subset(sub2) ~ sex + (1|p|fosternest) + (1|q|dam)) +
    set_rescor(FALSE)
  fit <- brm(bform, BTdata, refresh = 0)
  print(summary(fit))
  expect_error(predict(fit), "'resp' must be a single variable name")
  pred <- predict(fit, resp = "tarsus")
  expect_equal(nrow(pred), sum(BTdata$sub1))
  pred <- fitted(fit, resp = "back")
  expect_equal(nrow(pred), sum(BTdata$sub2))
  waic <- waic(fit, resp = "back")
  expect_range(waic$estimates[3, 1], 1100, 1200)
  ce <- conditional_effects(fit)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_equal(nobs(fit, resp = "tarsus"), sum(BTdata$sub1))
})

test_that("Cox models work correctly", {
  set.seed(12345)
  covs <- data.frame(id  = 1:200, trt = stats::rbinom(200, 1L, 0.5))
  d1 <- simsurv::simsurv(lambdas = 0.1, gammas  = 1.5, betas = c(trt = -0.5),
                         x = covs, maxt  = 5)
  d1 <- merge(d1, covs)

  fit1 <- brm(eventtime | cens(1 - status) ~ 1 + trt,
              data = d1, family = brmsfamily("cox"), refresh = 0)
  print(summary(fit1))
  expect_range(posterior_summary(fit1)["b_trt", "Estimate"], -0.70, -0.30)
  expect_range(waic(fit1)$estimates[3, 1], 620, 670)
})

test_that("ordinal model with grouped thresholds works correctly", {
  set.seed(1234)
  dat <- data.frame(
    y = sample(1:6, 100, TRUE),
    gr = rep(c("a", "b"), each = 50),
    th = rep(5:6, each = 50),
    x = rnorm(100)
  )

  prior <- prior(normal(0,1), class = "Intercept", group = "b")
  fit <- brm(y | thres(th, gr) ~ x, dat, cumulative(), prior = prior)
  print(summary(fit))
  pred <- predict(fit)
  expect_equal(dim(pred), c(nrow(dat), max(dat$th) + 1))
  expect_range(waic(fit)$estimates[3, 1], 350, 400)
  ce <- conditional_effects(fit, categorical = TRUE)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])

  # test incl_thres = TRUE
  thres_minus_eta <- posterior_linpred(fit, incl_thres = TRUE)
  bprep <- prepare_predictions(fit)
  thres <- bprep$thres$thres
  eta <- posterior_linpred(fit)
  gr_unq <- unique(family(fit)$thres$group)
  gr_vec <- fit$data$gr
  nthres_max <- max(
    by(family(fit)$thres, family(fit)$thres$group, function(x) max(x$thres))
  )
  thres_minus_eta_ch <- lapply(setNames(nm = gr_unq), function(gr) {
    thres_gr_nms <- grepl(paste0("^b_Intercept\\[", gr, ","), colnames(thres))
    thres_gr <- thres[, thres_gr_nms]
    eta_gr <- eta[, gr_vec == gr, drop = FALSE]
    thres_minus_eta_ch_gr <- apply(thres_gr, 2, "-", eta_gr)
    thres_minus_eta_ch_gr <- array(
      thres_minus_eta_ch_gr,
      dim = c(nrow(thres_gr), ncol(eta_gr), ncol(thres_gr))
    )
    if (ncol(thres_gr) < nthres_max) {
      dim_NA <- c(
        dim(thres_minus_eta_ch_gr)[-3],
        nthres_max - dim(thres_minus_eta_ch_gr)[3]
      )
      thres_minus_eta_ch_gr <-
        abind::abind(thres_minus_eta_ch_gr, array(dim = dim_NA))
    }
    dimnames(thres_minus_eta_ch_gr) <-
      list(NULL, NULL, as.character(seq_len(nthres_max)))
    return(thres_minus_eta_ch_gr)
  })
  new_arrnms <- dimnames(thres_minus_eta_ch[[1]])
  thres_minus_eta_ch <- abind::abind(thres_minus_eta_ch, along = 2)
  dimnames(thres_minus_eta_ch) <- new_arrnms
  expect_equivalent(thres_minus_eta, thres_minus_eta_ch)
})

test_that("Fixing parameters to constants works correctly", {
  bprior <- prior(normal(0, 1), class = "b") +
    prior(constant(2), class = "b", coef = "zBase") +
    prior(constant(0.5), class = "sd")

  fit <- brm(count ~ zAge + zBase + (1 | patient),
             data = epilepsy, prior = bprior)
  print(summary(fit))
  expect_range(waic(fit)$estimates[3, 1], 1790, 1840)
  ce <- conditional_effects(fit)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
})

test_that("projpred methods can be run", {
  fit <- brm(count ~ zAge + zBase * Trt,
             data = epilepsy, family = poisson())
  summary(fit)

  library(projpred)

  # perform variable selection without cross-validation
  vs <- varsel(fit)
  expect_is(vs, "vsel")

  # perform variable selection with cross-validation
  cv_vs <- cv_varsel(fit)
  expect_is(vs, "vsel")
})

test_that("emmeans can be run for multivariate models", {
  library(emmeans)
  df <- data.frame(
    y1 = rnorm(100), y2 = rnorm(100),
    x1 = rnorm(100), x2 = rnorm(100)
  )

  bform <- bf(mvbind(y1, y2) ~ x1 + x2) + set_rescor(TRUE)
  fit <- brm(bform, df, chains = 1, iter = 1000)

  # Default: Collapse over repeated measures:
  em <- summary(emmeans(fit, "x1", at = list(x1 = c(-1, 1))))
  expect_equal(nrow(em), 2)

  # Ask for MV with rep.meas
  em <- summary(emmeans(fit, c("x1", "rep.meas"), at = list(x1 = c(-1, 1))))
  expect_equal(nrow(em), 4)
})

test_that(paste(
  "Families sratio() and cratio() are equivalent for symmetric distribution",
  "functions (here only testing the logit link)"
), {
  set.seed(1234)
  dat2 <- data.frame(
    rating = sample(1:4, 50, TRUE),
    subject = rep(1:10, 5),
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50)
  )
  warmup <- 150
  iter <- 200
  chains <- 1

  fit_sratio <- brm(
    bf(rating ~ x1 + cs(x2) + (cs(x2)||subject), disc ~ 1),
    data = dat2, family = sratio(),
    warmup = warmup, iter = iter, chains = chains,
    seed = 533273
  )
  draws_sratio <- as.matrix(fit_sratio)

  fit_cratio <- brm(
    bf(rating ~ x1 + cs(x2) + (cs(x2)||subject), disc ~ 1),
    data = dat2, family = cratio(),
    warmup = warmup, iter = iter, chains = chains,
    seed = 533273
  )
  draws_cratio <- as.matrix(fit_cratio)

  expect_equal(draws_sratio, draws_cratio)
})

test_that("Non-linear non-looped model predictions work correctly in blocked order", {
  loss_alt <- transform(loss, row=as.integer(1:nrow(loss)), nr=nrow(loss), test=as.integer(0))
  scode_growth <- "
    vector growth_test(vector ult, int[] dev, vector theta, vector omega, int[] row, int[] test) {
      int N = rows(ult);
      vector[N] mu;
      int rows_sorted = 1;

      for(i in 1:N) {
        if(row[i] != i)
           rows_sorted = 0;
      }

      for(i in 1:N) {
         if(test[i] == 0) {
           mu[i] = ult[i] * (1 - exp(-(dev[i]/theta[i])^omega[i]));
         } else if(test[i] == 1) {
           mu[i] = rows_sorted;
         } else if(test[i] == 2) {
           mu[i] = N;
         }
      }
      return(mu);
    }
  "
  growth_model  <- stanvar(name="growth_test", scode=scode_growth, block="functions")

  fit_loss <- brm(
    bf(cum ~ growth_test(ult, dev, theta, omega, row, test),
       ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
       nl = TRUE, loop=FALSE),
    data = loss_alt, family = gaussian(),
    stanvars = growth_model,
    prior = c(prior(normal(5000, 1000), nlpar = "ult"),
              prior(normal(1, 2), nlpar = "omega"),
              prior(normal(45, 10), nlpar = "theta")),
    control = list(adapt_delta = 0.9),
    refresh = 0,
    chains = 1,
    backend = "rstan"
  )
  expose_functions(fit_loss)

  N <- nrow(loss_alt)
  pr1 <- posterior_epred(fit_loss, newdata=transform(loss_alt, test=1), ndraws=10)
  expect_true(all(dim(pr1) == c(10,N)))
  expect_true(all(pr1 == 1))
  pr1b <- posterior_epred(fit_loss, newdata=transform(loss_alt, test=1)[-5,], ndraws=10)
  expect_true(all(dim(pr1b) == c(10,N-1)))
  expect_true(all(pr1b == 0))
  pr1c <- posterior_epred(fit_loss, newdata=transform(loss_alt, test=1)[-N,], ndraws=10)
  expect_true(all(dim(pr1c) == c(10,N-1)))
  expect_true(all(pr1c == 1))
  pr2 <- posterior_epred(fit_loss, newdata=transform(loss_alt, test=2), ndraws=10)
  expect_true(all(dim(pr2) == c(10,N)))
  expect_true(all(pr2 == N))
  pr2b <- posterior_epred(fit_loss, newdata=transform(loss_alt, test=2)[-N,], ndraws=10)
  expect_true(all(dim(pr2b) == c(10,N-1)))
  expect_true(all(pr2b == N-1))
})

