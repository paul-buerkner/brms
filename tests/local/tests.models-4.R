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

test_that("Multiple imputation works correctly", suppressWarnings({
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
}))

test_that("Missing value imputation within Stan works correctly", suppressWarnings({
  data("nhanes", package = "mice")

  # add some new variables
  N <- nrow(nhanes)
  set.seed(5324)
  nhanes$sdy <- 5
  nhanes$sub <- TRUE
  nhanes$sub[1:2] <- FALSE
  nhanes$id <- 1:N
  nhanes$id[is.na(nhanes$chl)] <- c(500, 501)
  nhanes$idx <- sample(nhanes$id[nhanes$sub], N, TRUE)

  # basic missing value imputation
  bform1 <- bf(bmi | mi() ~ age * mi(chl)) +
    bf(chl | mi() ~ age) + set_rescor(FALSE)
  fit1 <- brm(bform1, data = nhanes, backend = "rstan", refresh = 0)

  print(fit1)
  pred1 <- predict(fit1)
  expect_true(!anyNA(pred1))
  ce1 <- conditional_effects(fit1, resp = "bmi")
  expect_ggplot(plot(ce1, ask = FALSE)[[1]])
  loo1 <- loo(fit1, newdata = na.omit(fit1$data))
  expect_range(loo1$estimates[3, 1], 200, 220)

  # missing value imputation with indexes
  bform2 <- bf(bmi | mi() ~ age * mi(chl, idx = idx)) +
    bf(chl | mi(idx = id) + subset(sub) ~ age) + set_rescor(FALSE)
  fit2 <- brm(bform2, data = nhanes, backend = "rstan", refresh = 0,
              control = list(adapt_delta = 0.99),
              iter = 3000, warmup = 1000)

  print(fit2)

  expect_true(all(c("Ymi_chl[500]", "Ymi_chl[501]") %in% variables(fit2)))
  expect_true(!any(paste0("Ymi_chl[", 1:25, "]") %in% variables(fit2)))

  pred2 <- predict(fit2, resp = "bmi")
  expect_true(!anyNA(pred2))

  ce2 <- conditional_effects(fit2, resp = "bmi")
  expect_ggplot(plot(ce2, ask = FALSE)[[1]])

  newdata <- nhanes
  newdata$bmi[is.na(newdata$bmi)] <- mean(newdata$bmi, na.rm = TRUE)
  loo2 <- loo(fit2, newdata = newdata, resp = "bmi")
  expect_range(loo2$estimates[3, 1], 260, 320)

  # overimputation
  bform3 <- bf(bmi | mi() ~ age * mi(chl)) +
    bf(chl | mi(sdy) ~ age) + set_rescor(FALSE)
  fit3 <- brm(bform3, data = nhanes,
              save_pars = save_pars(latent = TRUE),
              backend = "rstan", refresh = 0,
              control = list(adapt_delta = 0.99),
              iter = 3000, warmup = 1000)

  print(fit3)
  pred3 <- predict(fit3)
  expect_true(!anyNA(pred3))
  ce3 <- conditional_effects(fit3, resp = "bmi")
  expect_ggplot(plot(ce3, ask = FALSE)[[1]])
  loo3 <- loo(fit3, newdata = na.omit(fit3$data))
  expect_range(loo3$estimates[3, 1], 200, 230)

  # overimputation with indexes
  bform4 <- bf(bmi | mi() ~ age * mi(chl, idx = idx)) +
    bf(chl | mi(sdy, idx = id) + subset(sub) ~ age) + set_rescor(FALSE)
  fit4 <- brm(bform4, data = nhanes,
              save_pars = save_pars(latent = TRUE),
              backend = "rstan", refresh = 0,
              control = list(adapt_delta = 0.99),
              iter = 3000, warmup = 1000)

  print(fit4)

  id_values <- unique(nhanes$id[nhanes$sub])
  non_id_values <- unique(setdiff(nhanes$id[!nhanes$sub], id_values))
  expect_true(all(paste0("Ymi_chl[", id_values, "]") %in% variables(fit4)))
  expect_true(!any(paste0("Ymi_chl[", non_id_values, "]") %in% variables(fit4)))

  pred4 <- predict(fit4, resp = "bmi")
  expect_true(!anyNA(pred4))

  ce4 <- conditional_effects(fit4, resp = "bmi")
  expect_ggplot(plot(ce4, ask = FALSE)[[1]])

  newdata <- nhanes
  newdata$bmi[is.na(newdata$bmi)] <- mean(newdata$bmi, na.rm = TRUE)
  loo4 <- loo(fit4, newdata = newdata, resp = "bmi")
  expect_range(loo4$estimates[3, 1], 300, 350)
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
