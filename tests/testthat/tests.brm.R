test_that("Test that brm produces expected errors", {
  expect_error(brm(rating~treat+period+carry+(1|subject), data = inhaler, 
                   partial = ~treat, family = c("cratio", "logit")), 
              paste("Variables cannot be modeled as fixed", 
                    "and partial effects at the same time.", 
                    "Error occured for variables: treat"))
})

test_that("Test that all S3 methods have reasonable ouputs", {
  fit <- brmsfit_example
  # test S3 methods in alphabetical order
  # family
  expect_equal(family(fit), family("poisson", link = "log"))
  # fitted
  fitted_old <- fitted(fit)
  expect_equal(dim(fitted_old), c(nrow(epilepsy), 4))
  expect_equal(colnames(fitted_old), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  newdata <- data.frame(log_Age_c = c(0, -0.2), visit = c(1, 4),
                        Trt_c = c(-0.2, 0.5))
  fitted_new <- fitted(fit, newdata = newdata)
  expect_equal(dim(fitted_new), c(2, 4))
  # fixef
  fixef <- fixef(fit, estimate = c("mean", "sd"))  
  expect_equal(dimnames(fixef), list(c("Intercept", "log_Age_c"),
                                     c("mean", "sd")))
  # formula
  expect_equal(formula(fit), count ~ log_Age_c + (1|visit))
  # hypothesis
  h1 <- hypothesis(fit, "Intercept > log_Age_c")
  expect_equal(dim(h1$hypothesis), c(1, 6))
  expect_silent(print(h1))
  h2 <- hypothesis(fit, "Intercept = 0", class = "sd", group = "visit")
  expect_true(is.numeric(h2$hypothesis$Evid.Ratio[1]))
  expect_silent(print(h2))
  # omit launch_shiny
  # logLik
  expect_equal(dim(logLik(fit)), c(80, 236))
  # LOO
  .loo <- suppressWarnings(LOO(fit, cores = 1))
  expect_true(is.numeric(.loo[["looic"]]))
  expect_true(.loo[["se_looic"]] > 0)
  expect_silent(print(.loo))
  loo_compare2 <- suppressWarnings(LOO(fit, fit, cores = 1))
  expect_equal(length(loo_compare2), 2)
  expect_equal(dim(attr(loo_compare2, "compare")), c(1, 2))
  expect_silent(print(loo_compare2))
  loo_compare3 <- suppressWarnings(LOO(fit, fit, fit, cores = 1))
  expect_equal(length(loo_compare3), 3)
  expect_equal(dim(attr(loo_compare3, "compare")), c(3, 2))
  expect_silent(print(loo_compare3))
  # ngrps
  expect_equal(ngrps(fit), list(visit = 4))
  # nobs
  expect_equal(nobs(fit), nrow(epilepsy))
  # parnames 
  expect_equal(parnames(fit)[c(1, 3, 4, 10, 11)],
               c("b_Intercept", "sd_visit_Intercept", "r_visit[1]",
                 "prior_b", "lp__"))
  # plot tested in tests.plots.R
  # posterior_samples
  ps <- posterior_samples(fit)
  expect_equal(dim(ps), c(80, length(parnames(fit))))
  expect_equal(names(ps), parnames(fit))
  expect_equal(names(posterior_samples(fit, pars = "^b_")),
               c("b_Intercept", "b_log_Age_c"))
  # predict
  predict_old <- predict(fit)
  expect_equal(dim(predict_old), c(nrow(epilepsy), 4))
  expect_equal(colnames(predict_old), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  newdata <- data.frame(log_Age_c = c(0, -0.2), visit = c(1, 4))
  predict_new <- predict(fit, newdata = newdata)
  expect_equal(dim(predict_new), c(2, 4))
  # print
  expect_silent(print(fit))
  # prior_samples
  prs1 <- prior_samples(fit)
  expect_equal(dimnames(prs1),
               list(as.character(1:80),  c("sd_visit", "b_Intercept", "b")))
  prs2 <- prior_samples(fit, pars = "b_log_Age_c")
  expect_equal(dimnames(prs2), list(as.character(1:80), "b_log_Age_c"))
  expect_equal(prs1$b, prs2$b_log_Age_c)
  # ranef
  .ranef <- ranef(fit, estimate = "median", var = TRUE)
  expect_equal(dim(.ranef$visit), c(4, 1))
  expect_equal(dim(attr(.ranef$visit, "var")), c(1, 1, 4))
  # residuals
  res <- residuals(fit, type = "pearson", probs = c(0.65))
  expect_equal(dim(res), c(236, 3))
  # stancode
  expect_true(is.character(stancode(fit)))
  expect_silent(print(stancode(fit)))
  # standata
  expect_equal(names(standata(fit)),
               c("N", "Y", "K", "X", "J_1", "N_1", "K_1", "Z_1", "NC_1"))
  # stanplot tested in tests.plots.R
  # summary
  .summary <- summary(fit)
  expect_true(is.numeric(.summary$fixed))
  expect_equal(rownames(.summary$fixed), c("Intercept", "log_Age_c"))
  expect_equal(colnames(.summary$fixed), 
               c("Estimate", "Est.Error", "l-95% CI", 
                 "u-95% CI", "Eff.Sample", "Rhat"))
  expect_equal(rownames(.summary$random$visit), c("sd(Intercept)"))
  expect_true(is.numeric(.summary$WAIC))
  expect_silent(print(.summary))
  # do not test update as is causes CRAN checks to fail on Windows
  # VarCorr
  vc <- VarCorr(fit)
  expect_equal(names(vc), "visit")
  Names <- c("Intercept", "Trt_c")
  expect_equivalent(dimnames(vc$visit$cov$mean), 
                    list(Names, Names))
  expect_silent(print(vc))
  dat_vc <- as.data.frame(vc)
  expect_equal(dim(dat_vc), c(2, 7))
  expect_equal(names(data_vc), c("Estimate", "Group", "Name", "Std.Dev",
                                 "Cor", "Cov", "Cov"))
  # vcov
  expect_equal(dim(vcov(fit)), c(2, 2))
  # WAIC
  .waic <- WAIC(fit)
  expect_true(is.numeric(.waic[["waic"]]))
  expect_true(is.numeric(.waic[["se_waic"]]))
  waic_compare <- WAIC(fit, fit)
  expect_equal(length(waic_compare), 2)
  expect_equal(dim(attr(waic_compare, "compare")), c(1,2))
})