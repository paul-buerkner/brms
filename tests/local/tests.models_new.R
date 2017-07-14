rm(list = ls())
source("tests/local/setup.R")

test_that("Poisson model from brm doc works correctly", {
  ## Poisson regression for the number of seizures in epileptic patients
  ## using student_t priors for population-level effects
  ## and half cauchy priors for standard deviations of group-level effects
  fit1 <- brm(count ~ log_Age_c + log_Base4_c * Trt_c
              + (1|patient) + (1|obs),
              data = epilepsy, family = poisson(),
              prior = c(prior(student_t(5,0,10), class = b),
                        prior(cauchy(0,2), class = sd)),
              cores = 2)
  print(fit1)
  ## generate a summary of the results
  expect_range(fixef(fit1)[4, 1], -0.35, -0.30)
  ## extract random effects standard devations and covariance matrices
  dVC <- as.data.frame(VarCorr(fit1, estimate = c("mean", "sd"), old = TRUE))
  expect_equal(dim(dVC), c(4, 5))
  ## extract group specific effects of each level
  expect_equal(names(ranef(fit1)), c("obs", "patient"))
  ## predict responses based on the fitted model
  expect_equal(dim(predict(fit1)), c(nobs(fit1), 4))
  ## plot marginal effects of each predictor
  me1 <- marginal_effects(fit1)
  expect_ggplot(plot(me1, ask = FALSE)[[4]])
  ## investigare model fit
  expect_range(WAIC(fit1)$waic, 1120, 1160)
  expect_ggplot(pp_check(fit1))
  # test kfold
  kfold1 <- kfold(fit1, update_args = list(chains = 2, iter = 1000))
  expect_range(kfold1$kfoldic, 1220, 1260)
})

test_that("Ordinal model from brm doc works correctly", {
  ## Ordinal regression modeling patient's rating of inhaler instructions
  ## category specific effects are estimated for variable 'treat'
  fit2 <- brm(rating ~ period + carry + cs(treat),
              data = inhaler, family = sratio("cloglog"),
              prior = set_prior("normal(0,5)"),
              iter = 1000, chains = 2, cores = 2)
  print(fit2)
  expect_range(WAIC(fit2)$waic, 900, 950)
  expect_warning(me <- marginal_effects(fit2, effect = "treat"),
                 "Predictions are treated as continuous variables")
  expect_ggplot(plot(me)[[1]])
})

test_that("Survival model from brm doc works correctly", {
  ## Survival regression modeling the time between the first
  ## and second recurrence of an infection in kidney patients.
  fit3 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
              data = kidney, family = lognormal(),
              cores = 2)
  print(fit3)
  me3 <- marginal_effects(fit3, method = "predict")
  expect_ggplot(plot(me3, ask = FALSE)[[2]])
  expect_range(LOO(fit3)$looic, 650, 740)
})

test_that("Binomial model from brm doc works correctly", {
  ## Probit regression using the binomial family
  n <- sample(1:10, 100, TRUE)  # number of trials
  success <- rbinom(100, size = n, prob = 0.4)
  x <- rnorm(100)
  data4 <- data.frame(n, success, x)
  fit4 <- brm(success | trials(n) ~ x, data = data4,
              family = binomial("probit"),
              cores = 2)
  print(fit4)
  expect_message(me4 <- marginal_effects(fit4), 
                 "median number of trials")
  expect_ggplot(plot(me4, ask = FALSE)[[1]])
})

test_that("Models from hypotheses doc work correctly", {
  prior <- c(set_prior("normal(0,2)", class = "b"),
             set_prior("student_t(10,0,1)", class = "sigma"),
             set_prior("student_t(10,0,1)", class = "sd"))
  
  ## fit a linear mixed effects models
  fit <- brm(time ~ age + sex + disease + (1 + age|patient),
             data = kidney, family = lognormal(),
             prior = prior, sample_prior = TRUE,
             control = list(adapt_delta = 0.95),
             cores = 2)
  print(fit)
  
  ## perform two-sided hypothesis testing
  hyp1 <- hypothesis(fit, "sexfemale = age + diseasePKD")
  expect_range(hyp1$hypothesis$Estimate, 0.6, 0.74)
  expect_range(hyp1$hypothesis$Evid.Ratio, 2, 3.5)
  expect_true(is(plot(hyp1)[[1]], "ggplot"))
  
  ## perform one-sided hypothesis testing
  hyp2 <- hypothesis(fit, "diseasePKD + diseaseGN - 3 < 0")
  expect_range(hyp2$hypothesis$Evid.Ratio, 400)
  
  ## test more than one hypothesis at once
  hyp3 <- c("diseaseGN = diseaseAN", "2 * diseaseGN - diseasePKD = 0")
  hyp3 <- hypothesis(fit, hyp3)
  expect_equal(dim(hyp3$hypothesis), c(2, 6))
})

test_that("varying slopes without a fixed effect work", {
  fit1 <- brm(count ~ log_Age_c + log_Base4_c * Trt_c +
                (as.numeric(visit)|patient),
              data = epilepsy, family = gaussian(),
              chains = 2, cores = 2)
  print(fit1)
  me <- marginal_effects(fit1)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  # test reloo
  loo1 <- LOO(fit1)
  reloo1 <- reloo(loo1, fit1, chains = 1, iter = 100)
  expect_range(reloo1$looic, 1620, 1700)
  
  conditions <- data.frame(log_Age_c = 0, log_Base4_c = 0, Trt_c = 0)
  me <- marginal_effects(fit1, conditions = conditions)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  
  conditions <- data.frame(log_Age_c = 0, log_Base4_c = 0,
                           Trt_c = 0, visit = c(1:4, NA))
  me <- marginal_effects(fit1, conditions = conditions, re_formula = NULL)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  
  expect_range(WAIC(fit1)$waic, 1500, 1600)
  expect_equal(dim(predict(fit1)), c(nobs(fit1), 4))
})

test_that("categorical models work correctly", {
  fit2 <- brm(rating ~ period + carry + treat + (1|test|subject),
              data = inhaler, family = categorical, iter = 500,
              prior = c(set_prior("normal(0,5)", nlpar = "X2"),
                        set_prior("normal(0,5)", nlpar = "X3"),
                        set_prior("normal(0,5)", nlpar = "X4")),
              chains = 2, cores = 2)
  print(fit2)
  expect_range(WAIC(fit2)$waic, 830, 900)
  ncat <- length(unique(inhaler$rating))
  expect_equal(dim(predict(fit2)), c(nobs(fit2), ncat))
  expect_equal(dim(fitted(fit2)), c(nobs(fit2), 4, ncat))
})

test_that("ARMA models work correctly", {
  N <- 100
  y <- arima.sim(list(ar = c(0.7, -0.5, 0.04, 0.2, -0.4)), N)
  dat <- list(y = y, x = rnorm(N), g = sample(1:5, N, TRUE))
  
  fit_ar <- brm(y ~ x, data = dat, autocor = cor_ar(p = 5),
                prior = prior(normal(0, 5), class = "ar"),
                chains = 2, cores = 2)
  print(fit_ar)
  ar <- colMeans(as.matrix(fit_ar, "^ar"))
  expect_range(ar[1], 0.5, 0.9)
  expect_range(ar[3], -0.2, 0.25)
  expect_range(ar[5], -0.6, -0.1)
  expect_ggplot(plot(marginal_effects(fit_ar))[[1]])
  
  fit_ma <- brm(y ~ x, data = dat, autocor = cor_ma(q = 1),
                chains = 2, cores = 2)
  print(fit_ma)
  expect_gt(LOO(fit_ma)$looic, LOO(fit_ar)$looic)
  
  fit_arma <- brm(y ~ x + (1|g), data = dat,
                  autocor=cor_arma(~1|g, p = 1, q = 1, cov = TRUE),
                  prior = c(prior(normal(0, 5), class = "ar"),
                            prior(normal(0, 6), class = "ma")),
                  chains = 2, cores = 2)
  print(fit_arma)
  expect_range(waic(fit_arma)$waic, 280, 400)
  expect_equal(dim(predict(fit_arma)), c(nobs(fit_arma), 4))
  expect_ggplot(plot(marginal_effects(fit_arma))[[1]])
  
  fit_arr <- brm(y ~ x, data = dat, autocor = cor_arr(r = 5),
                 prior = prior(normal(0, 5), class = "arr"),
                 chains = 2, cores = 2)
})

test_that("multivariate normal models work correctly", {
  N <- 300
  y1 <- rnorm(N)
  y2 <- rnorm(N, mean = 1, sd = 2)
  x <- rnorm(N, 1)
  month <- sample(1:36, N, replace = TRUE)
  id <- sample(1:10, N, replace = TRUE)
  tim <- sample(1:N, 100)
  data <- data.frame(y1, y2, x, month, id, tim)
  
  fit_mv1 <- brm(cbind(y1,y2) ~ s(x) + poly(month, 3) + (1|x|id),
                  data = data, autocor = cor_arma(~tim|id, p = 1),
                  prior = c(prior_(~normal(0,5), resp = "y1"),
                            prior_(~normal(0,5), resp = "y2"),
                            prior_(~lkj(5), class = "rescor")),
                  sample_prior = TRUE,
                  iter = 1000, chains = 2, cores = 2)
  print(fit_mv1)
  
  dVC <- as.data.frame(VarCorr(fit_mv1, old = TRUE))
  expect_equal(dim(dVC), c(4, 7))
  expect_equal(dim(predict(fit_mv1)), c(600, 4))
  expect_equal(dim(fitted(fit_mv1)), c(600, 4))
  newdata <- data.frame(month = 1, y1 = 0, y2 = 0,
                        x = 0, id = 1, tim = 1)
  expect_equal(dim(predict(fit_mv1, newdata = newdata)), c(2, 4))
  ms <- marginal_smooths(fit_mv1, nsamples = 750)
  expect_equal(length(ms), 2)
  expect_ggplot(plot(ms, ask = FALSE)[[2]])
  
  fit_mv2 <- brm(cbind(y1, y2) ~ 1, data = data,
                 prior = prior_(~lkj(5), class = "rescor"),
                 sample_prior = TRUE, iter = 1000, cores = 2)
  print(fit_mv2)
  waic_mv <- WAIC(fit_mv1, fit_mv2, nsamples = 100)
  expect_equal(dim(waic_mv$ic_diffs__), c(1, 2))
})

test_that("ZI and HU models work correctly", {
  fit_hu <- brm(bf(count ~ log_Age_c + log_Base4_c * Trt_c + (1|id1|patient),
                   hu ~ log_Age_c + log_Base4_c * Trt_c + (1|id1|patient)),
                data = epilepsy, family = hurdle_poisson(),
                prior = c(prior(normal(0, 5)),
                          prior(normal(0, 5), nlpar = "hu")),
                chains = 2, cores = 2)
  print(fit_hu)
  expect_equal(dim(predict(fit_hu)), c(nobs(fit_hu), 4))
  expect_ggplot(plot(marginal_effects(fit_hu), ask = FALSE)[[2]])
  
  fit_zi <- brm(bf(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient),
                   zi ~ log_Age_c + log_Base4_c * Trt_c),
                data = epilepsy, family = zero_inflated_negbinomial(),
                prior = c(prior(normal(0,5)),
                          prior(normal(0,3), class = "sd"),
                          prior(cauchy(0,5), class = "shape")),
                chains = 2, cores = 2)
  print(fit_zi)
  expect_equal(dim(predict(fit_zi)), c(nobs(fit_zi), 4))
  expect_ggplot(plot(marginal_effects(fit_zi), ask = FALSE)[[2]])
  waic_zi <- WAIC(fit_hu, fit_zi, nsamples = 100)
  expect_equal(dim(waic_zi$ic_diffs__), c(1, 2))
  
  
  ## zero_inflated beta model
  data("GasolineYield", package = "betareg")
  dat <- GasolineYield
  dat$yield[c(1, 5, 8, 12, 16)] <- 0
  fit_zibeta <- brm(yield ~ batch + temp, data = dat,
                    family = zero_inflated_beta(),
                    chains = 2, cores = 2)
  print(fit_zibeta)
  expect_equal(dim(predict(fit_zibeta)), c(nobs(fit_zibeta), 4))
  expect_ggplot(plot(marginal_effects(fit_zibeta), ask = FALSE)[[1]])
  expect_range(WAIC(fit_zibeta)$waic, -100, -70)
})

test_that("Non-linear models work correctly", {
  load("tests/local/data/loss.Rda")
  head(loss)
  fit_loss <- brm(bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
                     ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
                     nl = TRUE),
                  data = loss, family = gaussian(),
                  prior = c(prior(normal(5000, 1000), nlpar = "ult"),
                            prior(normal(1, 2), nlpar = "omega"),
                            prior(normal(45, 10), nlpar = "theta")),
                  control = list(adapt_delta = 0.9))
  print(fit_loss)
  expect_ggplot(plot(marginal_effects(fit_loss))[[1]])
  expect_range(LOO(fit_loss)$looic, 700, 720)
})

test_that("Multivariate GAMMs work correctly", {
  n <- 200
  sig <- 2
  dat <- mgcv::gamSim(1, n = n, scale = sig)
  fit_gam <- brm(y ~ t2(x0, x2) + s(x1), data = dat,
                 chains = 2, cores = 2,
                 control = list(adapt_delta = 0.95))
  print(fit_gam)
  
  me <- marginal_effects(fit_gam)
  expect_ggplot(plot(me, ask = FALSE, rug = TRUE, points = TRUE)[[1]])
  me <- marginal_effects(fit_gam, surface = TRUE, too_far = 0.05)
  expect_ggplot(plot(me, ask = FALSE, rug = TRUE)[[1]])
  
  ms <- marginal_smooths(fit_gam, resolution = 25)
  expect_ggplot(plot(ms, rug = TRUE, ask = FALSE)[[1]])
  ms <- marginal_smooths(fit_gam, resolution = 100, too_far = 0.05)
  expect_ggplot(plot(ms, rug = TRUE, ask = FALSE)[[1]])
  
  expect_range(loo(fit_gam)$looic, 880, 940)
  expect_equal(dim(predict(fit_gam)), c(nobs(fit_gam), 4))
  
  newd <- data.frame(x0=(0:30)/30, x1=(0:30)/30,
                     x2=(0:30)/30, x3=(0:30)/30)
  prfi <- cbind(predict(fit_gam, newd), fitted(fit_gam, newdata = newd))
  expect_range(prfi[, 1], prfi[, 5] - 0.15, prfi[, 5] + 0.15)
})

test_that("GAMMs with factor variable in 'by' work correctly", {
  dat <- mgcv::gamSim(4, n = 200, dist = "normal")
  fit_gam2 <- brm(y ~ fac + s(x2, by = fac, k = 4), dat,
                  chains = 2, cores = 2)
  print(fit_gam2)
  
  me <- marginal_effects(fit_gam2, "x2:fac")
  expect_ggplot(plot(me, points = TRUE, ask = FALSE)[[1]])
  ms <- marginal_smooths(fit_gam2, res = 10)
  expect_ggplot(plot(ms, rug = TRUE, ask = FALSE)[[1]])
  
  fit_gam3 <- brm(y ~ fac + t2(x1, x2, by = fac), dat,
                  chains = 2, cores = 2)
  print(fit_gam3)
  
  me <- marginal_effects(fit_gam3, "x2:fac")
  expect_ggplot(plot(me, points = TRUE, ask = FALSE)[[1]])
  ms <- marginal_smooths(fit_gam3, too_far = 0.1)
  expect_ggplot(plot(ms, rug = TRUE, ask = FALSE)[[1]])
  expect_ggplot(plot(ms, rug = TRUE, stype = "raster", ask = FALSE)[[1]])
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
  
  fit_gev <- brm(bf(SeaLevel ~ cYear + SOI,
                    sigma ~ s(cYear, bs = "bs", m = 1, k = 3) + SOI),
                 data = fremantle,
                 family = gen_extreme_value(),
                 knots = knots, init_r = 0.5,
                 chains = 4, cores = 2,
                 control = list(adapt_delta = 0.95))
  print(fit_gev)
  
  prfi <- cbind(predict(fit_gev), fitted(fit_gev))
  expect_range(prfi[, 1], prfi[, 5] - 0.03, prfi[, 5] + 0.03)
  # expect_range(loo(fit_gev)$looic, -115, -95)
  me <- marginal_effects(fit_gev, "cYear")
  expect_ggplot(plot(me, points = TRUE, ask = FALSE)[[1]])
})

test_that("models are recompiled when changing number of FEs from 0 to 1", {
  fit1 <- brm(count ~ 1, data = epilepsy, cores = 2)
  fit2 <- update(fit1, ~ . + Trt_c, newdata = epilepsy, cores = 2)
  expect_equal(rownames(fixef(fit2)), c("Intercept", "Trt_c"))
  fit3 <- update(fit2, ~ . - Trt_c, newdata = epilepsy, cores = 2)
  expect_equal(rownames(fixef(fit3)), c("Intercept"))
})


test_that("family is correctly updated", {
  fit4 <- update(fit1, family = student())
  expect_equal(family(fit4)$family, "student")
  expect_equal(formula(fit4)$family$family, "student")
  fit5 <- update(fit1, bf(~., family = student()))
  expect_equal(family(fit5)$family, "student")
  expect_equal(formula(fit5)$family$family, "student")
})

test_that("Wiener diffusion models work correctly", {
  x <- rnorm(100, mean = 1)
  dat <- brms:::rwiener(n=1, alpha=2, tau=.3, beta=.5, delta=.5+x)
  dat$x <- x
  
  fit_d1 <- brm(bf(q | dec(resp) ~ x), dat, 
                family = wiener(), cores = 2)
  print(fit_d1)
  expect_ggplot(plot(marginal_effects(fit_d1), ask = FALSE)[[1]])
  expect_ggplot(pp_check(fit_d1))
  pp <- posterior_predict(fit_d1, nsamples = 10, negative_rt = TRUE)
  expect_true(min(pp) < 0)
  
  fit_d2 <- brm(bf(q | dec(resp) ~ x, ndt ~ x), 
                dat, family = wiener(), cores = 2)
  print(fit_d2)
  expect_ggplot(plot(marginal_effects(fit_d2), ask = FALSE)[[1]])
  expect_ggplot(pp_check(fit_d2))
  expect_warning(predict(fit_d2, nsamples = 10, negative_rt = TRUE),
                 "Summarizing positive and negative response times")
  
  waic_d <- WAIC(fit_d1, fit_d2)
  expect_equal(dim(waic_d$ic_diffs__), c(1, 2))
})

test_that("disc parameter in ordinal models is handled correctly", {
  fit <- brm(bf(rating ~ period + carry + treat + (1|subject), disc ~ 1),
             data = inhaler, family = cumulative(), 
             prior = c(prior(normal(0,5)),
                       prior(normal(0,1), "Intercept", nlpar = "disc")), 
             chains = 2, cores = 2)
  print(fit)
  expect_range(waic(fit)$waic, 870, 920)
  ncat <- length(unique(inhaler$rating))
  expect_equal(dim(predict(fit)), c(nobs(fit), ncat))
  expect_ggplot(plot(marginal_effects(fit), ask = FALSE, 
                     points = TRUE, jitter_width = 0.3)[[3]])
})

test_that("Mixture models work correctly", {
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
    prior(normal(0, 5), Intercept, nlpar = mu1),
    prior(normal(5, 5), Intercept, nlpar = mu2),
    prior(normal(10, 5), Intercept, nlpar = mu3),
    prior(dirichlet(1, 1, 1), theta)
  )
  fit1 <- brm(bform1, data = dat, family = mixfam, 
              prior = prior, chains = 2)
  print(fit1)
  expect_ggplot(pp_check(fit1))
  loo1 <- LOO(fit1)
  expect_equal(dim(pp_mixture(fit1)), c(nobs(fit1), 4, 3))
  
  bform2 <- bf(bform1, theta1 = 1, theta2 = 1, theta3 = 1)
  fit2 <- brm(bform2, data = dat, family = mixfam, 
              prior = prior, chains = 2)
  print(fit2)
  expect_ggplot(pp_check(fit2))
  loo2 <- LOO(fit2)
  expect_gt(loo2$looic, loo1$looic)
  expect_equal(dim(pp_mixture(fit2)), c(nobs(fit2), 4, 3))
  
  bform3 <- bf(bform1, theta1 ~ z, theta2 ~ 1)
  fit3 <- brm(bform3, data = dat, family = mixfam, 
              prior = prior, inits = 0, chains = 1)
  print(fit3)
  expect_ggplot(pp_check(fit3))
  loo3 <- LOO(fit3, pointwise = TRUE)
  expect_range(loo3$looic, loo1$looic - 20, loo1$looic + 20)
  expect_equal(dim(pp_mixture(fit3)), c(nobs(fit3), 4, 3))
})

test_that("Basic GPs work correctly", {
  ## Basic Gaussian process
  dat <- mgcv::gamSim(1, n = 30, scale = 2)
  fit1 <- brm(y ~ gp(x0) + x1 + gp(x2) + x3, dat, 
              chains = 2, cores = 2)
  print(fit1)
  expect_ggplot(pp_check(fit1))
  me = marginal_effects(fit1, nsamples = 200, nug = 1e-07)
  expect_ggplot(plot(me, ask = FALSE)[[3]])
  expect_range(WAIC(fit1)$waic, 100, 200)
})

test_that("Multivariate GPs work correctly", {
  fit2 <- brm(y ~ gp(x1, x2), dat, chains = 2, cores = 2)
  print(fit2)
  expect_ggplot(pp_check(fit2))
  me = marginal_effects(
    fit2, nsamples = 200, nug = 1e-07,
    surface = TRUE, resolution = 10
  )
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  expect_range(WAIC(fit2)$waic, 100, 200)
})
 
test_that("GPs with 'by' variable work correctly", {
  # continuous 'by' variable
  fit3 <- brm(y ~ gp(x1, by = x2), dat, chains = 2)
  print(fit3)
  expect_ggplot(pp_check(fit3))
  me = marginal_effects(fit3, nsamples = 200, nug = 1e-07)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  expect_range(WAIC(fit3)$waic, 100, 200)
  
  # factor 'by' variable
  dat2 <- mgcv::gamSim(4, n = 100, scale = 2)
  fit4 <- brm(y ~ gp(x2, by = fac), dat2, chains = 2, cores = 2)
  print(fit4)
  expect_ggplot(pp_check(fit4))
  me = marginal_effects(fit4, nsamples = 200, nug = 1e-07)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  expect_range(WAIC(fit4)$waic, 400, 600)
})

test_that("SAR models work correctly", {
  fit_lagsar <- brm(CRIME ~ INC + HOVAL, data= COL.OLD, 
                    family = student(),
                    autocor = cor_lagsar(COL.nb), 
                    chains = 2, cores = 2)
  print(fit_lagsar)
  expect_ggplot(pp_check(fit_lagsar))
  me = marginal_effects(fit_lagsar, nsamples = 200)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  expect_range(LOO(fit_lagsar)$loic, 350, 380)
  
  fit_errorsar <- brm(CRIME ~ INC + HOVAL, data = COL.OLD, 
                      autocor = cor_errorsar(COL.nb), 
                      chains = 2, cores = 2)
  print(fit_errorsar)
  expect_ggplot(pp_check(fit_errorsar))
  me = marginal_effects(fit_errorsar, nsamples = 200)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  expect_range(LOO(fit_errorsar)$looic, 350, 380)
})

test_that("CAR models work correctly", {
  # generate some spatial data
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
    y | trials(size) ~ x1 + x2, data = dat, 
    family = binomial(), autocor = cor_car(W, ~1|obs),
    chains = 2, cores = 2
  ) 
  print(fit_car)
  expect_ggplot(pp_check(fit_car))
  me = marginal_effects(fit_car, nsamples = 200)
  expect_ggplot(plot(me, ask = FALSE)[[1]])
  expect_range(LOO(fit_car)$looic, 450, 550)
  expect_false(isTRUE(all.equal(
    fitted(fit_car, newdata = dat[1:5, ]), 
    fitted(fit_car, newdata = dat[1:5, ], incl_autocor = FALSE)
  )))
  newdata <- data.frame(x1 = 0, x2 = 0, size = 50, obs = -1)
  expect_error(predict(fit_car, newdata = newdata),
               "Cannot handle new locations in CAR models")
})
