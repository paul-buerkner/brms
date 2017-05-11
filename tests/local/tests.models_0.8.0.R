rm(list = ls())

rerun <- FALSE
if (rerun) {
  packageurl <- "http://cran.r-project.org/src/contrib/Archive/brms/brms_0.8.0.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

source("tests/local/setup.R")

# old categorical models
if (rerun) {
  fit_old_cat <- brm(rating ~ period + carry + treat + (1|subject),
                     data = inhaler, family = categorical(),
                     prior = set_prior("normal(0,5)"),
                     chains = 2, iter = 500)
  # save fitted model objects
  save(fit_old_cat, file = "tests/local/models_0.8.0.Rda")
} else {
  load("tests/local/models_0.8.0.Rda")
}

print(fit_old_cat)
expect_range(WAIC(fit_old_cat)$waic, 870, 920)
expect_equal(dim(predict(fit_old_cat)), c(nobs(fit_old_cat), 4))
expect_equal(dim(fitted(fit_old_cat)), c(nobs(fit_old_cat), 4, 4))
