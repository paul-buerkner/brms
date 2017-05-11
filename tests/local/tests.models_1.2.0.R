rm(list = ls())

rerun <- FALSE
if (rerun) {
  packageurl <- "http://cran.r-project.org/src/contrib/Archive/brms/brms_1.2.0.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

source("tests/local/setup.R")


if (rerun) {
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
  # distributional model
  data_het <- data.frame(y = c(rnorm(50), rnorm(50, 1, 2)),
                         x = factor(rep(c("a", "b"), each = 50)))
  fit_old_dist <- brm(bf(y ~ x, sigma ~ 0 + x), data = data_het)
  # save fitted model objects
  save(fit_old_loss, fit_old_dist, file = "models_1.2.0.Rda")
} else {
  load("tests/local/models_1.2.0.Rda")
}

print(fit_old_loss)
expect_range(WAIC(fit_old_loss)$waic, 680, 730)
expect_equal(dim(predict(fit_old_loss)), c(nobs(fit_old_loss), 4))
expect_equal(dim(fitted(fit_old_loss)), c(nobs(fit_old_loss), 4))
expect_ggplot(plot(marginal_effects(fit_old_loss), ask = FALSE)[[1]])

print(fit_old_dist)
expect_range(WAIC(fit_old_dist)$waic, 350, 390)
expect_equal(dim(predict(fit_old_dist)), c(nobs(fit_old_dist), 4))
expect_equal(dim(fitted(fit_old_dist)), c(nobs(fit_old_dist), 4))
expect_ggplot(plot(marginal_effects(fit_old_dist), ask = FALSE)[[1]])
