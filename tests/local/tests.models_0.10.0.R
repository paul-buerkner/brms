rm(list = ls())

rerun <- FALSE
if (rerun) {
  packageurl <- "http://cran.r-project.org/src/contrib/Archive/brms/brms_0.10.0.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

source("tests/local/setup.R")

if (rerun) {
  # old MV linear model
  N <- 100
  y1 <- rnorm(N, sd = 5)
  y2 <- rnorm(N, sd = 2)
  y3 <- rnorm(N, 1)
  month <- sample(1:36, N, replace = TRUE)
  id <- sample(1:10, N, replace = TRUE)
  tim <- sample(1:100, N)
  data <- data.frame(y1 = y1, y2 = y2, y3 = y3, month = month, id = id, tim = tim)
  fit_old_mv<- brm(cbind(y1,y2,y3) ~ -1 + trait*poly(month,3) + (-1+trait|id), 
                   data = data, family = student(), 
                   autocor = cor_ar(~tim|id:trait), 
                   prior = c(set_prior("normal(0,5)"), 
                             set_prior("lkj(2)", class = "rescor")),
                   sample_prior = TRUE, chains = 2)
  # old ZI / HU model
  fit_old_hu <- brm(count ~ 0 + main + spec + 
                      (main + spec):(log_Age_c + log_Base4_c * Trt_c) +
                      (0+trait|visit),  
                    data = epilepsy, family = hurdle_poisson,
                    prior = set_prior("normal(0,5)"),
                    warmup = 500, iter = 1000, chains = 2)
  # old categorical model
  fit_old_cat <- brm(rating ~ 0 + trait + trait:(period + carry + treat) + 
                       (0+trait|subject), 
                     data = inhaler, family = categorical(), 
                     prior = set_prior("normal(0,5)"), 
                     chains = 2, iter = 500)
  # non-linear model
  url <- "https://raw.githubusercontent.com/mages/diesunddas/master/Data/ClarkTriangle.csv"
  loss <- read.csv(url)
  fit_old_loss <- brm(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
                      nonlinear = list(ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1), 
                      data = loss, family = gaussian(),
                      prior = c(set_prior("normal(5000, 1000)", nlpar = "ult"),
                                set_prior("normal(1, 2)", nlpar = "omega"),
                                set_prior("normal(45, 10)", nlpar = "theta")),
                      control = list(adapt_delta = 0.9))
  # save fitted model objects
  save(fit_old_mv, fit_old_hu, fit_old_cat, fit_old_loss, 
       file = "models_0.10.0.Rda")
} else {
  load("tests/local/models_0.10.0.Rda")
}

print(fit_old_mv)
expect_range(WAIC(fit_old_mv)$waic, 1320, 1360)
expect_equal(dim(predict(fit_old_mv)), c(nobs(fit_old_mv), 4))
expect_equal(dim(fitted(fit_old_mv)), c(nobs(fit_old_mv), 4))
expect_equal(dim(as.data.frame(VarCorr(fit_old_mv, old = TRUE))), c(6, 9))

print(fit_old_hu)
expect_range(WAIC(fit_old_hu)$waic, 1615, 1645)
expect_equal(dim(predict(fit_old_hu)), c(nobs(fit_old_hu) / 2, 4))
expect_equal(dim(fitted(fit_old_hu)), c(nobs(fit_old_hu) / 2, 4))
expect_ggplot(plot(marginal_effects(fit_old_hu), ask = FALSE)[[1]])

print(fit_old_cat)
expect_range(WAIC(fit_old_cat)$waic, 840, 865)
expect_equal(dim(predict(fit_old_cat)), c(nobs(fit_old_cat) / 3, 4))
expect_equal(dim(fitted(fit_old_cat)), c(nobs(fit_old_cat) / 3, 4, 4))

print(fit_old_loss)
expect_range(WAIC(fit_old_loss)$waic, 680, 730)
expect_equal(dim(predict(fit_old_loss)), c(nobs(fit_old_loss), 4))
expect_equal(dim(fitted(fit_old_loss)), c(nobs(fit_old_loss), 4))
expect_ggplot(plot(marginal_effects(fit_old_loss), ask = FALSE)[[1]])
