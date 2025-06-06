source("setup_tests_local.R")

test_that("ARMA models work correctly", suppressWarnings({
  set.seed(1234)
  N <- 100
  y <- arima.sim(list(ar = c(0.7, -0.5, 0.04, 0.2, -0.4)), N)
  dat <- list(y = y, x = rnorm(N), g = sample(1:5, N, TRUE))

  fit_ar <- brm(
    y ~ x + ar(p = 5),
    data = dat,
    prior = prior(normal(0, 5), class = "ar"),
    chains = 2, refresh = 0
  )
  print(fit_ar)
  ar <- colMeans(as.matrix(fit_ar, variable = "ar"))
  expect_range(ar[1], 0.5, 0.9)
  expect_range(ar[3], -0.2, 0.25)
  expect_range(ar[5], -0.6, -0.1)
  expect_ggplot(plot(conditional_effects(fit_ar))[[1]])

  fit_ma <- brm(y ~ x + ma(q = 1),
    data = dat,
    chains = 2, refresh = 0
  )
  print(fit_ma)
  expect_gt(LOO(fit_ma)$estimates[3, 1], LOO(fit_ar)$estimates[3, 1])

  fit_arma <- brm(
    y ~ x + (1 | g) + arma(gr = g, cov = TRUE),
    data = dat,
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
}))

test_that("categorical models work correctly", suppressWarnings({
  fit2 <- brm(rating ~ period + carry + treat + (1 | test | subject),
    data = inhaler, family = categorical, iter = 500,
    prior = c(
      set_prior("normal(0,5)", "b", dpar = c("mu2", "mu3", "mu4")),
      set_prior("normal(0,5)", "Intercept", dpar = c("mu2", "mu3", "mu4"))
    ),
    chains = 2, refresh = 0
  )
  print(fit2)
  expect_range(WAIC(fit2)$estimates[3, 1], 830, 900)
  ncat <- length(unique(inhaler$rating))
  expect_equal(dim(predict(fit2)), c(nobs(fit2), ncat))
  expect_equal(dim(fitted(fit2)), c(nobs(fit2), 4, ncat))
  expect_equal(
    dim(fitted(fit2, scale = "linear")),
    c(nobs(fit2), 4, ncat - 1)
  )

  # tests with new data
  newd <- inhaler[1:10, ]
  newd$rating <- NULL
  expect_equal(dim(predict(fit2, newdata = newd)), c(10, ncat))

  ce <- conditional_effects(fit2, categorical = TRUE)
  expect_ggplot(plot(ce, plot = FALSE)[[1]])
}))

test_that("multivariate normal models work correctly", suppressWarnings({
  set.seed(1234)
  N <- 300
  y1 <- rnorm(N)
  y2 <- rnorm(N, mean = 1, sd = 2)
  x <- rnorm(N, 1)
  month <- sample(1:36, N, replace = TRUE)
  id <- sample(1:10, N, replace = TRUE)
  tim <- sample(1:N, N)
  data <- data.frame(y1, y2, x, month, id, tim)

  fit_mv1 <- brm(
    mvbind(y1, y2) ~ s(x) + poly(month, 3) +
      (1 | x | id) + arma(tim, id, q = 0),
    data = data,
    prior = c(
      prior_(~ normal(0, 5), resp = "y1"),
      prior_(~ normal(0, 5), resp = "y2"),
      prior_(~ lkj(5), class = "rescor")
    ),
    sample_prior = TRUE,
    iter = 1000, chains = 2, refresh = 0
  )
  print(fit_mv1)

  expect_equal(dim(predict(fit_mv1)), c(300, 4, 2))
  expect_equal(dim(fitted(fit_mv1)), c(300, 4, 2))
  newdata <- data.frame(month = 1, y1 = 0, y2 = 0, x = 0, id = 1, tim = 1)
  expect_equal(dim(predict(fit_mv1, newdata = newdata)), c(1, 4, 2))
  cs <- conditional_smooths(fit_mv1, ndraws = 750)
  expect_equal(length(cs), 2)
  expect_ggplot(plot(cs, ask = FALSE)[[2]])

  fit_mv2 <- brm(mvbind(y1, y2) ~ 1,
    data = data,
    prior = prior_(~ lkj(5), class = "rescor"),
    sample_prior = TRUE, iter = 1000, refresh = 0
  )
  print(fit_mv2)
  waic_mv <- WAIC(fit_mv1, fit_mv2, ndraws = 100)
  expect_true(waic_mv$ic_diffs__[1, "WAIC"] > 0)
}))

test_that("emmeans can be run for multivariate models", suppressWarnings({
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
}))

test_that("generalized multivariate models work correctly", suppressWarnings({
  data("BTdata", package = "MCMCglmm")
  bform <- (bf(tarsus ~ sex + (1 | p | fosternest)) + skew_normal()) +
    (bf(back ~ s(tarsus, by = sex) + (1 | p | fosternest)) + gaussian())
  fit_mv <- brm(bform, BTdata, chains = 2, iter = 1000, refresh = 0)

  print(fit_mv)
  expect_ggplot(pp_check(fit_mv, resp = "back"))
  expect_range(waic(fit_mv)$estimates[3, 1], 4300, 4400)
  expect_ggplot(plot(conditional_effects(fit_mv), ask = FALSE)[[1]])
  expect_ggplot(plot(conditional_smooths(fit_mv))[[1]])
  expect_equal(dim(coef(fit_mv)$fosternest), c(104, 4, 7))
}))

test_that("ZI and HU models work correctly", suppressWarnings({
  fit_hu <- brm(
    bf(
      count ~ zAge + zBase * Trt + (1 | id1 | patient),
      hu ~ zAge + zBase * Trt + (1 | id1 | patient)
    ),
    data = epilepsy, family = hurdle_poisson(),
    prior = c(
      prior(normal(0, 5)),
      prior(normal(0, 5), dpar = "hu")
    ),
    chains = 2, refresh = 0
  )
  print(fit_hu)
  expect_equal(dim(predict(fit_hu)), c(nobs(fit_hu), 4))
  expect_ggplot(plot(conditional_effects(fit_hu), ask = FALSE)[[2]])

  fit_zi <- brm(
    bf(count ~ zAge + zBase * Trt + (1 | patient), zi ~ Trt),
    data = epilepsy, family = zero_inflated_negbinomial(),
    prior = prior(normal(0, 5)) +
      prior(normal(0, 3), class = "sd") +
      prior(cauchy(0, 5), class = "shape") +
      prior(normal(0, 5), dpar = "zi"),
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
    yield ~ batch + temp,
    data = dat,
    family = zero_inflated_beta(),
    chains = 2, inits = 0, refresh = 0
  )
  print(fit_zibeta)
  expect_equal(dim(predict(fit_zibeta)), c(nobs(fit_zibeta), 4))
  expect_ggplot(plot(conditional_effects(fit_zibeta), ask = FALSE)[[1]])
  expect_range(WAIC(fit_zibeta)$estimates[3, 1], -100, -70)
}))

test_that("Non-linear models work correctly", suppressWarnings({
  fit_loss <- brm(
    bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
      ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1,
      nl = TRUE
    ),
    data = loss, family = gaussian(),
    prior = c(
      prior(normal(5000, 1000), nlpar = "ult"),
      prior(normal(1, 2), nlpar = "omega"),
      prior(normal(45, 10), nlpar = "theta")
    ),
    control = list(adapt_delta = 0.9),
    refresh = 0
  )
  print(fit_loss)
  expect_ggplot(plot(conditional_effects(fit_loss))[[1]])
  expect_range(LOO(fit_loss)$estimates[3, 1], 700, 720)
}))

test_that(
  "Non-linear models of distributional parameters work correctly",
  suppressWarnings({
    set.seed(1234)
    x <- rnorm(100)
    y <- rnorm(100, 1, exp(x))
    dat <- data.frame(y, x, g = rep(1:10, each = 10))
    bform <- bf(y ~ x + (1 | V | g)) +
      nlf(sigma ~ a) +
      lf(a ~ x + (1 | V | g)) +
      gaussian()
    bprior <- prior(normal(0, 3), nlpar = a)
    fit <- brm(bform, dat, prior = bprior, chains = 2, refresh = 0)
    print(fit)
    ce <- conditional_effects(fit, method = "predict")
    expect_ggplot(plot(ce, ask = FALSE)[[1]])
    expect_equal(dim(fitted(fit, dat[1:10, ])), c(10, 4))
    expect_range(LOO(fit)$estimates[3, 1], 240, 350)
  })
)

test_that("Nested non-linear models work correctly", suppressWarnings({
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
}))

test_that(
  "Non-linear non-looped model predictions work correctly in blocked order",
  suppressWarnings({
    loss_alt <- transform(loss, row = as.integer(1:nrow(loss)), nr = nrow(loss), test = as.integer(0))
    scode_growth <- "
    vector growth_test(vector ult, array[] int dev, vector theta, vector omega, array[] int row, array[] int test) {
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
    growth_model <- stanvar(name = "growth_test", scode = scode_growth, block = "functions")

    fit_loss <- brm(
      bf(cum ~ growth_test(ult, dev, theta, omega, row, test),
        ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1,
        nl = TRUE, loop = FALSE
      ),
      data = loss_alt, family = gaussian(),
      stanvars = growth_model,
      prior = c(
        prior(normal(5000, 1000), nlpar = "ult"),
        prior(normal(1, 2), nlpar = "omega"),
        prior(normal(45, 10), nlpar = "theta")
      ),
      control = list(adapt_delta = 0.9),
      refresh = 0,
      chains = 1,
      backend = "rstan"
    )
    expose_functions(fit_loss)

    N <- nrow(loss_alt)
    pr1 <- posterior_epred(fit_loss, newdata = transform(loss_alt, test = 1), ndraws = 10)
    expect_true(all(dim(pr1) == c(10, N)))
    expect_true(all(pr1 == 1))
    pr1b <- posterior_epred(fit_loss, newdata = transform(loss_alt, test = 1)[-5, ], ndraws = 10)
    expect_true(all(dim(pr1b) == c(10, N - 1)))
    expect_true(all(pr1b == 0))
    pr1c <- posterior_epred(fit_loss, newdata = transform(loss_alt, test = 1)[-N, ], ndraws = 10)
    expect_true(all(dim(pr1c) == c(10, N - 1)))
    expect_true(all(pr1c == 1))
    pr2 <- posterior_epred(fit_loss, newdata = transform(loss_alt, test = 2), ndraws = 10)
    expect_true(all(dim(pr2) == c(10, N)))
    expect_true(all(pr2 == N))
    pr2b <- posterior_epred(fit_loss, newdata = transform(loss_alt, test = 2)[-N, ], ndraws = 10)
    expect_true(all(dim(pr2b) == c(10, N - 1)))
    expect_true(all(pr2b == N - 1))
  })
)
