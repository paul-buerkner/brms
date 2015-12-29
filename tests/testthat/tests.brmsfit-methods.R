test_that("all S3 methods have reasonable ouputs", {
  fit <- rename_pars(brmsfit_example)
  # test S3 methods in alphabetical order
  # family
  expect_equal(family(fit), family("poisson", link = "log"))
  # fitted
  fitted1 <- fitted(fit)
  expect_equal(dim(fitted1), c(nrow(epilepsy), 4))
  expect_equal(colnames(fitted1), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  
  newdata <- data.frame(log_Age_c = c(0, -0.2), visit = c(1, 4),
                        Trt_c = c(-0.2, 0.5))
  
  fitted2 <- fitted(fit, newdata = newdata)
  expect_equal(dim(fitted2), c(2, 4))
  newdata$visit <- c(1, 6)
  
  fitted3 <- fitted(fit, newdata = newdata, 
                    allow_new_levels = TRUE)
  expect_equal(dim(fitted3), c(2, 4))
  # fixef
  fixef <- fixef(fit, estimate = c("mean", "sd"))  
  expect_equal(dimnames(fixef), list(c("Intercept", "Trt_c"),
                                     c("mean", "sd")))
  # formula
  expect_equal(formula(fit), 
               count ~ Trt_c + offset(log_Age_c) + (1+Trt_c|visit))
  # hypothesis
  h1 <- hypothesis(fit, "Intercept > Trt_c")
  expect_equal(dim(h1$hypothesis), c(1, 6))
  expect_silent(capture.output(print(h1)))
  expect_silent(p <- plot(h1, do_plot = FALSE))
  
  h2 <- hypothesis(fit, "Intercept = 0", class = "sd", group = "visit")
  expect_true(is.numeric(h2$hypothesis$Evid.Ratio[1]))
  expect_silent(capture.output(print(h2)))
  expect_silent(p <- plot(h2, ignore_prior = TRUE, do_plot = FALSE))
  expect_error(hypothesis(fit, "Intercept > x"), fixed = TRUE,
               "cannot be found in the model: b_x")
  # omit launch_shiny
  # logLik
  expect_equal(dim(logLik(fit)), c(80, 236))
  # LOO
  .loo <- suppressWarnings(LOO(fit, cores = 1))
  expect_true(is.numeric(.loo[["looic"]]))
  expect_true(.loo[["se_looic"]] > 0)
  expect_silent(capture.output(print(.loo)))
  
  loo_compare2 <- suppressWarnings(LOO(fit, fit, cores = 1))
  expect_equal(length(loo_compare2), 2)
  expect_equal(dim(attr(loo_compare2, "compare")), c(1, 2))
  expect_silent(capture.output(print(loo_compare2)))
  
  loo_compare3 <- suppressWarnings(LOO(fit, fit, fit, cores = 1))
  expect_equal(length(loo_compare3), 3)
  expect_equal(dim(attr(loo_compare3, "compare")), c(3, 2))
  expect_silent(capture.output(print(loo_compare3)))
  # model.frame
  expect_equal(model.frame(fit), fit$data)
  # ngrps
  expect_equal(ngrps(fit), list(visit = 4))
  # nobs
  expect_equal(nobs(fit), nrow(epilepsy))
  # parnames 
  expect_equal(parnames(fit)[c(1, 3, 7, 16, 18)],
               c("b_Intercept", "sd_visit_Intercept", "r_visit[2,1]",
                 "prior_b", "lp__"))
  # plot tested in tests.plots.R
  # posterior_samples
  ps <- posterior_samples(fit)
  expect_equal(dim(ps), c(80, length(parnames(fit))))
  expect_equal(names(ps), parnames(fit))
  expect_equal(names(posterior_samples(fit, pars = "^b_")),
               c("b_Intercept", "b_Trt_c"))
  # predict
  predict1 <- predict(fit)
  expect_equal(dim(predict1), c(nrow(epilepsy), 4))
  expect_equal(colnames(predict1), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  expect_equal(dim(predict(fit, nsamples = 10, probs = 0.5)), 
               c(nrow(epilepsy), 3))
  
  newdata <- data.frame(log_Age_c = c(0, -0.2), visit = c(1, 4),
                        Trt_c = c(-0.2, 0.5))
  predict2 <- predict(fit, newdata = newdata)
  expect_equal(dim(predict2), c(2, 4))
  
  newdata$visit <- c(1, 6)
  predict3 <- predict(fit, newdata = newdata, 
                      allow_new_levels = TRUE)
  expect_equal(dim(predict3), c(2, 4))
  # print
  expect_silent(capture.output(print(fit)))
  # prior_samples
  prs1 <- prior_samples(fit)
  expect_equal(dimnames(prs1),
               list(as.character(1:80), 
                    c("sd_visit", "b_Intercept", "b", "cor_visit")))
  
  prs2 <- prior_samples(fit, pars = "b_Trt_c")
  expect_equal(dimnames(prs2), list(as.character(1:80), "b_Trt_c"))
  expect_equal(prs1$b, prs2$b_Trt_c)
  # ranef
  .ranef <- ranef(fit, estimate = "median", var = TRUE)
  expect_equal(dim(.ranef$visit), c(4, 2))
  expect_equal(dim(attr(.ranef$visit, "var")), c(2, 2, 4))
  # residuals
  res1 <- residuals(fit, type = "pearson", probs = c(0.65))
  expect_equal(dim(res1), c(236, 3))
  newdata <- epilepsy[1:10, ]
  
  res2 <- residuals(fit, newdata = newdata)
  expect_equal(dim(res2), c(10, 4))
  newdata$visit <- rep(1:5, 2)
  
  res3 <- residuals(fit, newdata = newdata,
                    allow_new_levels = TRUE)
  expect_equal(dim(res3), c(10, 4))
  # stancode
  expect_true(is.character(stancode(fit)))
  expect_silent(capture.output(print(stancode(fit))))
  # standata
  expect_equal(names(standata(fit)),
               c("N", "Y", "offset", "K", "X", 
                 "J_1", "N_1", "K_1", "Z_1", "NC_1"))
  # stanplot tested in tests.plots.R
  # summary
  .summary <- summary(fit)
  expect_true(is.numeric(.summary$fixed))
  expect_equal(rownames(.summary$fixed), c("Intercept", "Trt_c"))
  expect_equal(colnames(.summary$fixed), 
               c("Estimate", "Est.Error", "l-95% CI", 
                 "u-95% CI", "Eff.Sample", "Rhat"))
  expect_equal(rownames(.summary$random$visit), 
               c("sd(Intercept)", "sd(Trt_c)", "cor(Intercept,Trt_c)"))
  expect_true(is.numeric(.summary$WAIC))
  expect_silent(capture.output(print(.summary)))
  # update
  # do not actually refit the model as is causes CRAN checks to fail
  new_data <- data.frame(log_Age_c = c(0, 1, -1), visit = c(3, 2, 4),
                         Trt_c = c(0, 0.5, -0.5))
  up <- update(fit, newdata = new_data, ranef = FALSE, refit = FALSE)
  expect_true(class(up) == "brmsfit")
  expect_equal(up$data.name, "new_data")
  expect_equal(attr(up$ranef$visit, "levels"), c("2", "3", "4"))
  expect_true("r_1" %in% up$exclude)
  expect_error(update(fit, family = "gaussian"),
               "family cannot be updated")
  expect_error(update(fit, data = new_data),
               "use argument 'newdata' to update your data")
  # VarCorr
  vc <- VarCorr(fit)
  expect_equal(names(vc), "visit")
  Names <- c("Intercept", "Trt_c")
  expect_equivalent(dimnames(vc$visit$cov$mean), 
                    list(Names, Names))
  expect_silent(capture.output(print(vc)))
  data_vc <- as.data.frame(vc)
  expect_equal(dim(data_vc), c(2, 7))
  expect_equal(names(data_vc), c("Estimate", "Group", "Name", 
                                 "Std.Dev", "Cor", "Cov", "Cov"))
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