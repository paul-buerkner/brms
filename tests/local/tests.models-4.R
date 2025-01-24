source("setup_tests_local.R")

test_that("UNSTR models work correctly", suppressWarnings({
  epilepsy2 <- epilepsy
  epilepsy2$visit <- as.numeric(epilepsy2$visit)
  fit <- brm(count ~ Trt + unstr(visit, patient), data = epilepsy2)
  print(fit)
  expect_ggplot(pp_check(fit))

  waic <- waic(fit)
  expect_range(waic$estimates[3, 1], 1550, 1600)
  # ensure the the correlation are actually included in the predictions
  waic_without <- waic(fit, incl_autocor = FALSE)
  expect_true(waic$estimates[3, 1] + 200 < waic_without$estimates[3, 1])

  waic_new <- waic(fit, newdata = epilepsy2[1:100, ])
  expect_range(waic_new $estimates[3, 1], 700, 760)

  newdat <- epilepsy2[1:5, ]
  newdat$visit[1] <- 5
  expect_error(
    waic(fit, newdata = newdat),
    "Cannot handle new time points in UNSTR models"
  )
}))

test_that("SAR models work correctly", suppressWarnings({
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
}))

test_that("CAR models work correctly", suppressWarnings({
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
}))

test_that("Missing value imputation works correctly", suppressWarnings({
  library(mice)
  data("nhanes", package = "mice")

  # missing value imputation via multiple imputation
  imp <- mice(nhanes, m = 5)
  fit_imp1 <- brm_multiple(bmi ~ age * chl, imp, chains = 1,
                           iter = 2000, warmup = 1000,
                           backend = "rstan", refresh = 0)
  print(fit_imp1)
  expect_equal(ndraws(fit_imp1), 5000)

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
  expect_range(loo$estimates[3, 1], 200, 225)
}))

test_that("student-t-distributed group-level effects work correctly",
          suppressWarnings({
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
}))

test_that("multinomial models work correctly", suppressWarnings({
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
}))

test_that("dirichlet_multinomial models work correctly", suppressWarnings({
  require("extraDistr")
  set.seed(1245)
  N <- 100
  dat <- as.data.frame(extraDistr::rdirmnom(N, 10, c(10, 5, 1)))
  names(dat) <- paste0("y", 1:3)
  dat$size <- with(dat, y1 + y2 + y3)
  dat$x <- rnorm(N)
  dat$y <- with(dat, cbind(y1, y2, y3))
  
  fit <- brm(
    y | trials(size) ~ x, data = dat,
    family = dirichlet_multinomial(),
    prior = prior("exponential(0.01)", "phi")
  )
  print(summary(fit))
  pred <- predict(fit)
  expect_equal(dim(pred), c(nobs(fit), 4, 3))
  expect_equal(dimnames(pred)[[3]], c("y1", "y2", "y3"))
  waic <- waic(fit)
  expect_range(waic$estimates[3, 1], 550, 650)
  ce <- conditional_effects(fit, categorical = TRUE)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
}))

test_that("dirichlet models work correctly", suppressWarnings({
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
}))

test_that("logistic_normal models work correctly", suppressWarnings({
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
}))
