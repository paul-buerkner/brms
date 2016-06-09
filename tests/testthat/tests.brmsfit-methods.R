test_that("all S3 methods have reasonable ouputs", {
  fit <- rename_pars(brmsfit_example)
  # test S3 methods in alphabetical order
  # as.data.frame
  ps <- as.data.frame(fit)
  expect_true(is(ps, "data.frame"))
  expect_equal(dim(ps), c(Nsamples(fit), length(parnames(fit))))
  # as.matrix
  ps <- as.matrix(fit)
  expect_true(is(ps, "matrix"))
  expect_equal(dim(ps), c(Nsamples(fit), length(parnames(fit))))
  # as.mcmc
  chains <- fit$fit@sim$chains
  mc <- as.mcmc(fit)
  expect_equal(length(mc), chains)
  expect_equal(dim(mc[[1]]), c(Nsamples(fit) / chains, length(parnames(fit))))
  # test assumes thin = 1
  expect_equal(dim(as.mcmc(fit, inc_warmup = TRUE)[[1]]), 
               c(fit$fit@sim$iter, length(parnames(fit))))
  # coef
  expect_equal(dim(coef(fit)$visit), c(4, 3))
  # family
  expect_equal(family(fit), family("gaussian", link = "identity"))
  # fitted
  fitted1 <- fitted(fit)
  expect_equal(dim(fitted1), c(nrow(epilepsy), 4))
  expect_equal(colnames(fitted1), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  
  newdata <- data.frame(log_Age_c = c(0, -0.2), visit = c(1, 4),
                        Trt_c = c(-0.2, 0.5), count = c(20, 13),
                        patient = c(1, 42), Exp = c(2, 4))
  fitted2 <- fitted(fit, newdata = newdata)
  expect_equal(dim(fitted2), c(2, 4))
  newdata$visit <- c(1, 6)
  
  fitted3 <- fitted(fit, newdata = newdata, 
                    allow_new_levels = TRUE)
  expect_equal(dim(fitted3), c(2, 4))
  # fixef
  fixef <- fixef(fit, estimate = c("mean", "sd"))  
  expect_equal(dimnames(fixef), list(c("Intercept", "Trt_c", "Exp"),
                                     c("mean", "sd")))
  # formula
  expect_equal(formula(fit), 
     count ~ Trt_c + mono(Exp) + offset(log_Age_c) + (1+Trt_c|visit))
  # hypothesis
  h1 <- hypothesis(fit, "Intercept > Trt_c")
  expect_equal(dim(h1$hypothesis), c(1, 6))
  expect_output(print(h1), "Intercept-(Trt_c) > 0", fixed = TRUE)
  expect_silent(p <- plot(h1, do_plot = FALSE))
  
  h2 <- hypothesis(fit, "Intercept = 0", class = "sd", group = "visit")
  expect_true(is.numeric(h2$hypothesis$Evid.Ratio[1]))
  expect_output(print(h2), "class sd_visit:", fixed = TRUE)
  expect_silent(p <- plot(h2, ignore_prior = TRUE, do_plot = FALSE))
  expect_error(hypothesis(fit, "Intercept > x"), fixed = TRUE,
               "cannot be found in the model: b_x")
  # omit launch_shiny
  # logLik
  expect_equal(dim(logLik(fit)), c(Nsamples(fit), nobs(fit)))
  # LOO
  .loo <- suppressWarnings(LOO(fit, cores = 1))
  expect_true(is.numeric(.loo[["looic"]]))
  expect_true(.loo[["se_looic"]] > 0)
  expect_output(print(.loo), "LOOIC")
  
  loo_compare2 <- suppressWarnings(LOO(fit, fit, cores = 1))
  expect_equal(length(loo_compare2), 2)
  expect_equal(dim(attr(loo_compare2, "compare")), c(1, 2))
  expect_output(print(loo_compare2), "fit - fit")
  
  loo_compare3 <- suppressWarnings(LOO(fit, fit, fit, cores = 1))
  expect_equal(length(loo_compare3), 3)
  expect_equal(dim(attr(loo_compare3, "compare")), c(3, 2))
  #expect_output(print(loo_compare3), "Weights")
  # marginal_effects (the related plot method is tested in tests.plots)
  mdata = data.frame(log_Age_c = c(-0.3, 0, 0.3), count = c(10, 20, 30), 
                     visit = 1:3, patient = 1, Trt_c = 0, Exp = c(1,3,5))
  exp_nrow <- nrow(mdata) * 100
  expect_equal(nrow(marginal_effects(fit, conditions = mdata)[[1]]),
               exp_nrow)
  expect_equal(nrow(marginal_effects(fit, effects = "Trt_c", 
                                     conditions = mdata)[[1]]), 
               exp_nrow)
  expect_equal(nrow(marginal_effects(fit, re_formula = NULL, 
                                     conditions = mdata)[[1]]), 
               exp_nrow)
  expect_error(marginal_effects(fit), "Please specify argument 'conditions' manually")
  expect_error(marginal_effects(fit, effects = "Trt_cc"), 
               "All specified effects are invalid for this model")
  # model.frame
  expect_equal(model.frame(fit), fit$data)
  # ngrps
  expect_equal(ngrps(fit), list(visit = 4))
  # nobs
  expect_equal(nobs(fit), nrow(epilepsy))
  # parnames 
  expect_equal(parnames(fit)[c(1, 3, 4, 8, 11, 21, 23, 31)],
               c("b_Intercept", "b_Exp", "ar[1]", "cor_visit_Intercept_Trt_c", 
                 "simplex_Exp[2]", "r_visit[4,Trt_c]", "prior_sigma", "lp__"))
  # plot tested in tests.plots.R
  # posterior_samples
  ps <- posterior_samples(fit)
  expect_equal(dim(ps), c(Nsamples(fit), length(parnames(fit))))
  expect_equal(names(ps), parnames(fit))
  expect_equal(names(posterior_samples(fit, pars = "^b_")),
               c("b_Intercept", "b_Trt_c", "b_Exp"))
  # predict
  predict1 <- predict(fit)
  expect_equal(dim(predict1), c(nrow(epilepsy), 4))
  expect_equal(colnames(predict1), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  expect_equal(dim(predict(fit, nsamples = 10, probs = 0.5)), 
               c(nrow(epilepsy), 3))
  
  newdata <- data.frame(log_Age_c = c(0, -0.2), visit = c(1, 4),
                        Trt_c = c(-0.2, 0.5), count = c(2, 10),
                        patient = c(1, 42), Exp = c(1, 2))
  predict2 <- predict(fit, newdata = newdata)
  expect_equal(dim(predict2), c(2, 4))
  
  newdata$visit <- c(1, 6)
  predict3 <- predict(fit, newdata = newdata, 
                      allow_new_levels = TRUE)
  expect_equal(dim(predict3), c(2, 4))
  # print
  expect_output(SW(print(fit)), "Group-Level Effects:")
  # prior_samples
  prs1 <- prior_samples(fit)
  prior_names <- c("sd_visit", "sigma", "b", "bm", 
                   paste0("simplex_Exp[", 1:4, "]"), "cor_visit")
  expect_equal(dimnames(prs1),
               list(as.character(1:Nsamples(fit)), prior_names))
  
  prs2 <- prior_samples(fit, pars = "b_Trt_c")
  expect_equal(dimnames(prs2), list(as.character(1:Nsamples(fit)), "b_Trt_c"))
  expect_equal(sort(prs1$b), sort(prs2$b_Trt_c))
  # ranef
  .ranef <- ranef(fit, estimate = "median", var = TRUE)
  expect_equal(dim(.ranef$visit), c(4, 2))
  expect_equal(dim(attr(.ranef$visit, "var")), c(2, 2, 4))
  # residuals
  res1 <- residuals(fit, type = "pearson", probs = c(0.65))
  expect_equal(dim(res1), c(nobs(fit), 3))
  newdata <- cbind(epilepsy[1:10, ], Exp = rep(1:5, 2))
  res2 <- residuals(fit, newdata = newdata)
  expect_equal(dim(res2), c(10, 4))
  newdata$visit <- rep(1:5, 2)
  
  res3 <- residuals(fit, newdata = newdata,
                    allow_new_levels = TRUE)
  expect_equal(dim(res3), c(10, 4))
  # stancode
  expect_true(is.character(stancode(fit)))
  expect_output(print(stancode(fit)), "generated quantities")
  # standata
  expect_equal(names(standata(fit)),
               c("N", "Y", "K", "X_means", "X", "Km", "Xm", "Jm", 
                 "con_simplex_1", "J_1", "N_1", "K_1", "NC_1", 
                 "Z_1_1", "Z_1_2", "offset", "tg", "Kar", "Kma", 
                 "Karma", "prior_only"))
  # stanplot tested in tests.plots.R
  # summary
  .summary <- SW(summary(fit, waic = TRUE))
  expect_true(is.numeric(.summary$fixed))
  expect_equal(rownames(.summary$fixed), c("Intercept", "Trt_c", "Exp"))
  expect_equal(colnames(.summary$fixed), 
               c("Estimate", "Est.Error", "l-95% CI", 
                 "u-95% CI", "Eff.Sample", "Rhat"))
  expect_equal(rownames(.summary$random$visit), 
               c("sd(Intercept)", "sd(Trt_c)", "cor(Intercept,Trt_c)"))
  expect_true(is.numeric(.summary$WAIC))
  expect_output(print(.summary), "Population-Level Effects:")
  # update
  # do not actually refit the model as is causes CRAN checks to fail
  new_data <- data.frame(log_Age_c = c(0, 1, -1), visit = c(3, 2, 4),
                         Trt_c = c(0, 0.5, -0.5), count = c(5, 17, 28),
                         patient = 1, Exp = 4)
  up <- update(fit, newdata = new_data, ranef = FALSE, testmode = TRUE)
  expect_true(class(up) == "brmsfit")
  expect_equal(up$data.name, "new_data")
  expect_equal(attr(up$ranef$visit, "levels"), c("2", "3", "4"))
  expect_true("r_1" %in% up$exclude)
  expect_error(update(fit, data = new_data), "use argument 'newdata'")
  # VarCorr
  vc <- VarCorr(fit)
  expect_equal(names(vc), c("visit", "RESIDUAL"))
  Names <- c("Intercept", "Trt_c")
  expect_equivalent(dimnames(vc$visit$cov$mean), 
                    list(Names, Names))
  expect_output(print(vc), "visit")
  data_vc <- as.data.frame(vc)
  expect_equal(dim(data_vc), c(3, 7))
  expect_equal(names(data_vc), c("Estimate", "Group", "Name", 
                                 "Std.Dev", "Cor", "Cov", "Cov"))
  # vcov
  expect_equal(dim(vcov(fit)), c(3, 3))
  # WAIC
  .waic <- WAIC(fit)
  expect_true(is.numeric(.waic[["waic"]]))
  expect_true(is.numeric(.waic[["se_waic"]]))
  waic_compare <- WAIC(fit, fit)
  expect_equal(length(waic_compare), 2)
  expect_equal(dim(attr(waic_compare, "compare")), c(1,2))
})