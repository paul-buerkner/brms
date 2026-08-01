source("setup_tests_local.R")

test_that("Multivariate GAMMs work correctly", suppressWarnings({
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
  cs <- conditional_smooths(fit_gam, surface = FALSE)
  expect_ggplot(plot(cs, rug = TRUE, ask = FALSE)[[1]])

  expect_range(loo(fit_gam)$estimates[3, 1], 830, 870)
  expect_equal(dim(predict(fit_gam)), c(nobs(fit_gam), 4))

  newd <- data.frame(x0 = (0:30)/30, x1 = (0:30)/30,
                     x2 = (0:30)/30, x3 = (0:30)/30)
  prfi <- cbind(predict(fit_gam, newd), fitted(fit_gam, newdata = newd))
  expect_range(prfi[, 1], prfi[, 5] - 0.25, prfi[, 5] + 0.25)
}))

test_that("GAMMs with factor variable in 'by' work correctly", suppressWarnings({
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
}))

test_that("generalized extreme value models work correctly", suppressWarnings({
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
    knots = knots, init = 0.5, chains = 4,
    control = list(adapt_delta = 0.95), refresh = 0
  )
  print(fit_gev)

  prfi <- cbind(predict(fit_gev), fitted(fit_gev))
  expect_range(prfi[, 1], prfi[, 5] - 0.03, prfi[, 5] + 0.03)
  # expect_range(loo(fit_gev)$estimates[3, 1], -115, -95)
  ce <- conditional_effects(fit_gev, "cYear")
  expect_ggplot(plot(ce, points = TRUE, ask = FALSE)[[1]])
}))

test_that("Wiener diffusion models work correctly", suppressWarnings({
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
}))

test_that("disc parameter in ordinal models is handled correctly", suppressWarnings({
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
}))

test_that("Argument `incl_thres` works correctly for non-grouped thresholds",
          suppressWarnings({
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
}))

test_that("hurdle_cumulative family works correctly", suppressWarnings({
  inhaler2 <- inhaler
  inhaler2$rating[1:10] <- 0
  fit <- brm(
    bf(rating ~ period + carry + treat, hu ~ treat),
    data = inhaler2, family = hurdle_cumulative(),
    prior = prior(normal(0,5)),
    chains = 2, refresh = 0
  )
  print(fit)
  expect_range(waic(fit)$estimates[3, 1], 950, 1050)
  ncat <- length(unique(inhaler$rating))
  expect_equal(dim(predict(fit)), c(nobs(fit), ncat))
  expect_ggplot(plot(
    SW(conditional_effects(fit)), ask = FALSE,
    points = TRUE, point_args = list(width = 0.3)
  )[[3]])
}))

test_that("Mixture models work correctly", suppressWarnings({
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
    chains = 2, init = 0, refresh = 0, seed = 1234
  )
  print(fit1)
  expect_ggplot(pp_check(fit1))
  loo1 <- loo(fit1)
  expect_equal(dim(pp_mixture(fit1)), c(nobs(fit1), 4, 3))

  bform2 <- bf(bform1, theta1 = 1, theta2 = 1, theta3 = 1)
  fit2 <- brm(
    bform2, data = dat, family = mixfam,
    prior = prior, chains = 2,
    init = 0, refresh = 0
  )
  print(fit2)
  expect_ggplot(pp_check(fit2))
  loo2 <- loo(fit2)
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
}))

test_that("group-level mixture models work correctly", suppressWarnings({
  set.seed(1234)
  ng <- 40
  npg <- 6
  comp <- rbinom(ng, 1, 0.4)
  dat <- do.call(rbind, lapply(seq_len(ng), function(g) {
    data.frame(g = factor(g), y = rnorm(npg, c(0, 6)[comp[g] + 1], 1))
  }))
  fit <- brm(
    bf(y ~ 1), data = dat,
    family = mixture(gaussian, gaussian, gr = "g"),
    prior = c(prior(normal(0, 5), Intercept, dpar = mu1),
              prior(normal(6, 5), Intercept, dpar = mu2)),
    chains = 2, init = 0, refresh = 0, seed = 1234,
    save_pars = save_pars(all = TRUE)  # required for loo_moment_match
  )
  print(fit)
  # the mixing proportions and component means are recovered
  ps <- posterior_summary(fit, variable = c("b_mu1_Intercept", "b_mu2_Intercept"))
  expect_range(ps["b_mu1_Intercept", "Estimate"], -1, 1)
  expect_range(ps["b_mu2_Intercept", "Estimate"], 5, 7)
  # predictions are per observation, but log_lik / loo are per group
  expect_equal(dim(posterior_predict(fit)), c(ndraws(fit), nobs(fit)))
  expect_equal(dim(posterior_epred(fit)), c(ndraws(fit), nobs(fit)))
  expect_equal(ncol(log_lik(fit)), ng)
  expect_equal(nrow(loo(fit)$pointwise), ng)
  # pointwise evaluation is also computed per group
  expect_equal(
    loo(fit, pointwise = TRUE)$pointwise[, "elpd_loo"],
    loo(fit)$pointwise[, "elpd_loo"]
  )
  expect_equal(dim(pp_mixture(fit)), c(ng, 4, 2))
  # kfold and the loo refinements operate per group (leave-one-group-out)
  loo1 <- loo(fit)
  # exact leave-one-group-out kfold tracks PSIS-loo
  kf <- kfold(fit, folds = "loo", chains = 1, refresh = 0)
  expect_equal(nrow(kf$pointwise), ng)
  expect_equal(
    kf$estimates["elpd_kfold", "Estimate"],
    loo1$estimates["elpd_loo", "Estimate"], tolerance = 2
  )
  # K-fold over groups keeps whole groups within a fold
  kf5 <- kfold(fit, K = 5, chains = 1, refresh = 0, save_fits = TRUE)
  expect_equal(nrow(kf5$pointwise), ng)
  # predictions from the saved fits are per observation
  kfp <- kfold_predict(kf5)
  expect_equal(length(kfp$y), nobs(fit))
  # a single integer 'Ksub' selects that many random folds
  kf2 <- kfold(fit, K = 5, Ksub = 2, chains = 1, refresh = 0)
  expect_lt(nrow(kf2$pointwise), ng)
  # reloo refits problematic groups and preserves the per-group unit
  rl <- reloo(fit, loo1, k_threshold = 0.5, chains = 1, refresh = 0)
  expect_equal(nrow(rl$pointwise), ng)
  expect_true(all(is.finite(rl$pointwise[, "elpd_loo"])))
  # subsampling and moment matching operate on the per-group pointwise term
  expect_equal(nrow(loo_subsample(fit, observations = 20)$pointwise), 20)
  expect_s3_class(loo_moment_match(fit, loo1, k_threshold = 0.5), "loo")
  # 'group' is redundant with the mixture grouping and is rejected
  expect_error(kfold(fit, group = "g"), "not supported for group-level")
  # newdata may contain a subset of the original groups or new groups
  nd <- dat[dat$g %in% c("1", "2"), ]
  nd$g <- factor(rep(c("1", "new"), each = npg))
  expect_equal(ncol(log_lik(fit, newdata = nd)), 2)
  expect_equal(
    dim(posterior_predict(fit, newdata = nd)),
    c(ndraws(fit), nrow(nd))
  )
  expect_equal(dim(pp_mixture(fit, newdata = nd)), c(2, 4, 2))
}))

test_that("Gaussian processes work correctly", suppressWarnings({
  ## Basic GPs
  set.seed(1112)
  dat <- mgcv::gamSim(1, n = 30, scale = 2)
  fit1 <- brm(y ~ gp(x0) + x1 + gp(x2, cov = "matern32") + x3, dat,
              chains = 2, refresh = 0)
  print(fit1)
  expect_ggplot(pp_check(fit1))
  ce <- conditional_effects(fit1, ndraws = 200, nug = 1e-07)
  expect_ggplot(plot(ce, ask = FALSE)[[3]])
  expect_range(WAIC(fit1)$estimates[3, 1], 100, 200)

  # multivariate isotropic GPs
  fit2 <- brm(y ~ gp(x1, x2), dat, chains = 2, refresh = 0)
  print(fit2)
  expect_ggplot(pp_check(fit2))
  ce <- conditional_effects(
    fit2, ndraws = 200, nug = 1e-07,
    surface = TRUE, resolution = 10
  )
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit2)$estimates[3, 1], 100, 200)

  # multivariate non-isotropic GPs
  fit3 <- brm(y ~ gp(x1, x2, iso = FALSE, cov = "matern32"),
              data = dat, chains = 2, refresh = 0)
  print(fit3)
  expect_ggplot(pp_check(fit3))
  ce <- conditional_effects(
    fit3, ndraws = 200, nug = 1e-07,
    surface = TRUE, resolution = 10
  )
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit3)$estimates[3, 1], 70, 130)

  # GP with continuous 'by' variable
  fit4 <- brm(y ~ gp(x1, by = x2, cov = "matern52"), dat,
              chains = 2, refresh = 0)
  print(fit4)
  expect_ggplot(pp_check(fit4))
  ce <- conditional_effects(fit4, ndraws = 200, nug = 1e-07)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit4)$estimates[3, 1], 100, 200)

  # GP with factor 'by' variable
  dat2 <- mgcv::gamSim(4, n = 100, scale = 2)
  fit5 <- brm(y ~ gp(x2, by = fac, cov = "exponential"), dat2,
              chains = 2, refresh = 0)
  print(fit5)
  expect_ggplot(pp_check(fit5))
  ce <- conditional_effects(fit5, ndraws = 200, nug = 1e-07)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_range(WAIC(fit5)$estimates[3, 1], 400, 600)
}))

test_that("Approximate Gaussian processes work correctly", suppressWarnings({
  set.seed(1245)
  dat <- mgcv::gamSim(4, n = 200, scale = 2)

  # isotropic approximate GP
  fit1 <- brm(
    y ~ gp(x1, x2, by = fac, k = 10, c = 5/4, cov = "matern32"),
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

  # non-isotropic approximate GP
  fit2 <- brm(
    y ~ gp(x1, x2, by = fac, k = 10, c = 5/4, iso = FALSE,
           cov = "matern52"),
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
}))
