test_that("all S3 methods have reasonable ouputs", {
  fit1 <- brms:::rename_pars(brms:::brmsfit_example1)
  fit2 <- brms:::rename_pars(brms:::brmsfit_example2)
  fit3 <- brms:::rename_pars(brms:::brmsfit_example3)
  fit4 <- brms:::rename_pars(brms:::brmsfit_example4)
  
  # test S3 methods in alphabetical order
  # as.data.frame
  ps <- as.data.frame(fit1)
  expect_true(is(ps, "data.frame"))
  expect_equal(dim(ps), c(Nsamples(fit1), length(parnames(fit1))))
  
  # as.matrix
  ps <- as.matrix(fit1)
  expect_true(is(ps, "matrix"))
  expect_equal(dim(ps), c(Nsamples(fit1), length(parnames(fit1))))
  
  # as.mcmc
  chains <- fit1$fit@sim$chains
  mc <- as.mcmc(fit1)
  expect_equal(length(mc), chains)
  expect_equal(dim(mc[[1]]), c(Nsamples(fit1) / chains, length(parnames(fit1))))
  mc <- as.mcmc(fit1, combine_chains = TRUE)
  expect_equal(dim(mc), c(Nsamples(fit1), length(parnames(fit1))))
  # test assumes thin = 1
  expect_equal(dim(as.mcmc(fit1, inc_warmup = TRUE)[[1]]), 
               c(fit1$fit@sim$iter, length(parnames(fit1))))
  
  # coef
  coef1 <- coef(fit1)
  expect_equal(dim(coef1$visit), c(4, 8))
  coef2 <- coef(fit2)
  expect_equal(dim(coef2[[2]]), c(59, 2))
  expect_equal(attr(coef2[[1]], "nlpar"), "a")
  coef4 <- coef(fit4)
  expect_equal(dim(coef4$subject), c(10, 8))
  
  # family
  expect_equal(family(fit1), brmsfamily("student", link = "identity"))
  
  # fitted
  fi <- fitted(fit1)
  expect_equal(dim(fi), c(nobs(fit1), 4))
  expect_equal(colnames(fi), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  
  newdata <- data.frame(Age = c(0, -0.2), visit = c(1, 4),
                        Trt = c(-0.2, 0.5), count = c(20, 13),
                        patient = c(1, 42), Exp = c(2, 4))
  fi <- fitted(fit1, newdata = newdata)
  expect_equal(dim(fi), c(2, 4))
  newdata$visit <- c(1, 6)
  fi <- fitted(fit1, newdata = newdata, 
               allow_new_levels = TRUE)
  expect_equal(dim(fi), c(2, 4))

  fi <- fitted(fit2)
  expect_equal(dim(fi), c(nobs(fit2), 4))
  fi <- fitted(fit2, newdata = newdata,
               allow_new_levels = TRUE)
  expect_equal(dim(fi), c(2, 4))
  
  fi <- fitted(fit4)
  expect_equal(dim(fi), c(nobs(fit4), 4, 4))
  
  # fixef
  fixef1 <- fixef(fit1, estimate = c("mean", "sd"))  
  expect_equal(dimnames(fixef1), 
               list(c("Intercept", "Trt", "Age", "Trt:Age", "sAge_1", 
                      "sigma_Intercept", "sigma_Trt", "Exp"),
                    c("mean", "sd")))
  fixef2 <- fixef(fit2, estimate = c("median", "quantile"), 
                  probs = c(0.05, 0.95))
  expect_equal(dimnames(fixef2),
               list(c("a_Intercept", "a_Age", "b_Intercept", "b_Age"),
                    c("median", "5%", "95%")))
  
  # formula
  expect_equal(formula(fit1), 
    count ~ Trt * Age + mono(Exp) + s(Age) + offset(Age) + (1 + Trt | visit))
  
  # hypothesis
  h1 <- hypothesis(fit1, "Intercept > Trt")
  expect_equal(dim(h1$hypothesis), c(1, 6))
  expect_output(print(h1), "(Intercept)-(Trt) > 0", fixed = TRUE)
  expect_silent(p <- plot(h1, plot = FALSE))
  
  h2 <- hypothesis(fit1, "Intercept = 0", class = "sd", group = "visit")
  expect_true(is.numeric(h2$hypothesis$Evid.Ratio[1]))
  expect_output(print(h2), "class sd_visit:", fixed = TRUE)
  expect_silent(p <- plot(h2, ignore_prior = TRUE, plot = FALSE))
  expect_error(hypothesis(fit1, "Intercept > x"), fixed = TRUE,
               "cannot be found in the model: \nb_x")
  
  # omit launch_shiny
  
  # log_lik
  expect_equal(dim(log_lik(fit1)), c(Nsamples(fit1), nobs(fit1)))
  expect_equal(dim(log_lik(fit2)), c(Nsamples(fit2), nobs(fit2)))
  expect_equal(log_lik(fit1), logLik(fit1))
  
  # LOO
  loo1 <- SW(LOO(fit1, cores = 1))
  expect_true(is.numeric(loo1[["looic"]]))
  expect_true(loo1[["se_looic"]] > 0)
  expect_output(print(loo1), "LOOIC")
  expect_equal(loo1, SW(loo(fit1, cores = 1)))
  
  loo_compare1 <- SW(LOO(fit1, fit1, cores = 1))
  expect_equal(length(loo_compare1), 2)
  expect_equal(dim(attr(loo_compare1, "compare")), c(1, 2))
  expect_output(print(loo_compare1), "fit1 - fit1")
  
  loo_compare2 <- SW(LOO(fit1, fit1, fit1, cores = 1))
  expect_equal(length(loo_compare2), 3)
  expect_equal(dim(attr(loo_compare2, "compare")), c(3, 2))
  # expect_output(print(loo_compare3), "Weights")
  
  loo2 <- SW(LOO(fit2, cores = 1))
  expect_true(is.numeric(loo2[["looic"]]))
  
  loo3 <- SW(LOO(fit3, cores = 1))
  expect_true(is.numeric(loo3[["looic"]]))
  loo3 <- SW(LOO(fit3, pointwise = TRUE, cores = 1))
  expect_true(is.numeric(loo3[["looic"]]))
  
  loo4 <- SW(LOO(fit4, cores = 1))
  expect_true(is.numeric(loo4[["looic"]]))
  
  # marginal_effects (the related plot method is tested in tests.plots)
  expect_equal(nrow(marginal_effects(fit1)[[2]]), 100)
  mdata = data.frame(Age = c(-0.3, 0, 0.3), count = c(10, 20, 30), 
                     visit = 1:3, patient = 1, Trt = 0, Exp = c(1,3,5))
  exp_nrow <- nrow(mdata) * 100
  expect_equal(nrow(marginal_effects(fit1, conditions = mdata)[[1]]),
               exp_nrow)
  expect_equal(nrow(marginal_effects(fit1, effects = "Trt", 
                                     conditions = mdata)[[1]]), 
               exp_nrow)
  expect_equal(nrow(marginal_effects(fit1, re_formula = NULL, 
                                     conditions = mdata)[[1]]), 
               exp_nrow)
  expect_error(marginal_effects(fit1, effects = "Trtc"), 
               "All specified effects are invalid for this model")
  expect_warning(marginal_effects(fit1, effects = c("Trtc", "Trt")), 
                 "Some specified effects are invalid for this model")
  expect_error(marginal_effects(fit1, effects = "Trtc:a:b"), 
               "please use the 'conditions' argument")
  expect_equal(nrow(marginal_effects(fit2)[[2]]), 100)
  expect_equal(nrow(marginal_effects(fit2, conditions = mdata)[[1]]),
               exp_nrow)
  expect_warning(me4 <- marginal_effects(fit4),
                 "Predictions are treated as continuous variables")
  expect_true(is(me4, "brmsMarginalEffects"))
  
  # marginal_smooths
  ms1 <- marginal_smooths(fit1)
  expect_equal(nrow(ms1[[1]]), 100)
  expect_true(is(ms1, "brmsMarginalEffects"))
  expect_error(marginal_smooths(fit1, smooths = "s3"),
               "No valid smooth terms found in the model")
  expect_error(marginal_smooths(fit2),
               "No valid smooth terms found in the model")
  
  # model.frame
  expect_equal(model.frame(fit1), fit1$data)
  
  # ngrps
  expect_equal(ngrps(fit1), list(visit = 4))
  expect_equal(ngrps(fit2), list(patient = 59))
  
  # nobs
  expect_equal(nobs(fit1), nrow(epilepsy))
  
  # parnames 
  expect_equal(parnames(fit1)[c(1, 8, 9, 13, 15, 17, 27, 35, 38, 46)],
               c("b_Intercept", "bmo_Exp", "ar[1]", "cor_visit__Intercept__Trt", 
                 "nu", "simplex_Exp[2]", "r_visit[4,Trt]", "s_sAge_1[8]", 
                 "prior_sd_visit", "lp__"))
  expect_equal(parnames(fit2)[c(1, 4, 6, 7, 9, 71, 129)],
               c("b_a_Intercept", "b_b_Age", "sd_patient__b_Intercept",
                 "cor_patient__a_Intercept__b_Intercept", 
                 "r_patient__a[1,Intercept]", "r_patient__b[4,Intercept]",
                 "prior_b_a"))
  
  # plot tested in tests.plots.R
  
  # posterior_samples
  ps <- posterior_samples(fit1)
  expect_equal(dim(ps), c(Nsamples(fit1), length(parnames(fit1))))
  expect_equal(names(ps), parnames(fit1))
  expect_equal(names(posterior_samples(fit1, pars = "^b_")),
               c("b_Intercept", "b_Trt", "b_Age", "b_Trt:Age", 
                 "b_sAge_1", "b_sigma_Intercept", "b_sigma_Trt"))
  
  # posterior_predict
  expect_equal(dim(posterior_predict(fit1)), 
               c(Nsamples(fit1), nobs(fit1)))
  
  # pp_check
  expect_true(is(pp_check(fit1), "ggplot"))
  expect_true(is(pp_check(fit1, newdata = fit1$data[1:100, ]), "ggplot"))
  expect_true(is(pp_check(fit1, "stat", nsamples = 5), "ggplot"))
  expect_true(is(pp_check(fit1, "error_binned"), "ggplot"))
  ribbon_plot <- pp_check(fit1, "ribbon_grouped", group = "visit", x = "Age")
  expect_true(is(ribbon_plot, "ggplot"))
  expect_true(is(pp_check(fit3), "ggplot"))
  expect_true(is(pp_check(fit2, "ribbon", x = "Trt"), "ggplot"))
  expect_error(pp_check(fit1, "wrong_type"))
  expect_error(pp_check(fit2, "violin_grouped"), "group")
  expect_error(pp_check(fit1, "stat_grouped", group = "g"),
               "not a valid grouping factor")
  expect_true(is(pp_check(fit4), "ggplot"))
  
  # predict
  pred <- predict(fit1)
  expect_equal(dim(pred), c(nobs(fit1), 4))
  expect_equal(colnames(pred), 
               c("Estimate", "Est.Error", "2.5%ile", "97.5%ile"))
  expect_equal(dim(predict(fit1, nsamples = 10, probs = 0.5)), 
               c(nobs(fit1), 3))
  
  newdata <- data.frame(Age = c(0, -0.2), visit = c(1, 4),
                        Trt = c(-0.2, 0.5), count = c(2, 10),
                        patient = c(1, 42), Exp = c(1, 2))
  pred <- predict(fit1, newdata = newdata)
  expect_equal(dim(pred), c(2, 4))
  
  newdata$visit <- c(1, 6)
  pred <- predict(fit1, newdata = newdata, 
                  allow_new_levels = TRUE)
  expect_equal(dim(pred), c(2, 4))
  
  pred <- predict(fit2)
  expect_equal(dim(pred), c(nobs(fit2), 4))
  
  pred <- predict(fit2, newdata = newdata, 
                  allow_new_levels = TRUE)
  expect_equal(dim(pred), c(2, 4))
  
  pred <- predict(fit4)
  expect_equal(dim(pred), c(nobs(fit4), 4))
  # check if grouping factors with a single level are accepted
  newdata$patient <- factor(2)
  pred <- predict(fit2, newdata = newdata)
  expect_equal(dim(pred), c(2, 4))
  
  # predictive error
  expect_equal(dim(predictive_error(fit1)), 
               c(Nsamples(fit1), nobs(fit1)))
  
  # print
  expect_output(SW(print(fit1)), "Group-Level Effects:")
  
  # prior_samples
  prs1 <- prior_samples(fit1)
  prior_names <- c("sds_sAge_1", "nu", "sd_visit", "b", "bmo", 
                   paste0("simplex_Exp[", 1:4, "]"), "cor_visit")
  expect_equal(dimnames(prs1),
               list(as.character(1:Nsamples(fit1)), prior_names))
  
  prs2 <- prior_samples(fit1, pars = "b_Trt")
  expect_equal(dimnames(prs2), list(as.character(1:Nsamples(fit1)), "b_Trt"))
  expect_equal(sort(prs1$b), sort(prs2$b_Trt))
  
  # prior_summary
  expect_true(is(prior_summary(fit1), "brmsprior"))
  
  # ranef
  ranef1 <- ranef(fit1, estimate = "median", var = TRUE)
  expect_equal(dim(ranef1$visit), c(4, 2))
  expect_equal(dim(attr(ranef1$visit, "var")), c(2, 2, 4))
  
  ranef2 <- ranef(fit2, estimate = "mean", var = TRUE)
  expect_equal(dim(ranef2[[1]]), c(59, 1))
  expect_equal(dim(attr(ranef2[[1]], "var")), c(1, 1, 59))
  
  # residuals
  res1 <- residuals(fit1, type = "pearson", probs = c(0.65))
  expect_equal(dim(res1), c(nobs(fit1), 3))
  newdata <- cbind(epilepsy[1:10, ], Exp = rep(1:5, 2))
  res2 <- residuals(fit1, newdata = newdata)
  expect_equal(dim(res2), c(10, 4))
  newdata$visit <- rep(1:5, 2)
  
  res3 <- residuals(fit1, newdata = newdata,
                    allow_new_levels = TRUE)
  expect_equal(dim(res3), c(10, 4))
  
  res4 <- residuals(fit2)
  expect_equal(dim(res4), c(nobs(fit2), 4))
  
  expect_error(residuals(fit4), 
               "Residuals not implemented for family 'sratio'")
  
  # stancode
  expect_true(is.character(stancode(fit1)))
  expect_output(print(stancode(fit1)), "generated quantities")
  
  # standata
  expect_equal(names(standata(fit1)),
               c("N", "Y", "nb_1", "knots_1", "Zs_1_1", "K", "X", 
                 "Kmo", "Xmo", "Jmo", "con_simplex_1", "Z_1_1", "Z_1_2", 
                 "offset", "K_sigma", "X_sigma", "J_1", "N_1", "M_1", 
                 "NC_1", "tg", "Kar", "Kma", "Karma", "prior_only"))
  expect_equal(names(standata(fit2)),
               c("N", "Y", "KC", "C", "K_a", "X_a", "Z_1_a_1",
                 "K_b", "X_b", "Z_1_b_2", "J_1", "N_1", "M_1",
                 "NC_1", "disp", "prior_only"))
  
  # stanplot tested in tests.plots.R
  # summary
  summary1 <- SW(summary(fit1, waic = TRUE, priors = TRUE))
  expect_true(is.numeric(summary1$fixed))
  expect_equal(rownames(summary1$fixed), 
               c("Intercept", "Trt", "Age", "Trt:Age", "sAge_1", 
                 "sigma_Intercept", "sigma_Trt", "Exp"))
  expect_equal(colnames(summary1$fixed), 
               c("Estimate", "Est.Error", "l-95% CI", 
                 "u-95% CI", "Eff.Sample", "Rhat"))
  expect_equal(rownames(summary1$random$visit), 
               c("sd(Intercept)", "sd(Trt)", "cor(Intercept,Trt)"))
  expect_true(is.numeric(summary1$WAIC))
  expect_output(print(summary1), "Population-Level Effects:")
  expect_output(print(summary1), "Priors:")
  
  summary2 <- SW(summary(fit1, waic = TRUE))
  
  # update
  # do not actually refit the model as is causes CRAN checks to fail
  up <- update(fit1, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  new_data <- data.frame(Age = rnorm(18), visit = rep(c(3, 2, 4), 6),
                         Trt = rep(c(0, 0.5, -0.5), 6), 
                         count = rep(c(5, 17, 28), 6),
                         patient = 1, Exp = 4)
  up <- update(fit1, newdata = new_data, ranef = FALSE, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  expect_equal(up$data.name, "new_data")
  expect_equal(attr(up$ranef, "levels")$visit, c("2", "3", "4"))
  expect_true("r_1" %in% up$exclude)
  expect_error(update(fit1, data = new_data), "use argument 'newdata'")
  up <- update(fit1, formula = ~ . + log(Trt), testmode = TRUE,
               prior = set_prior("normal(0,10)"))
  expect_true(is(up, "brmsfit"))
  up <- update(fit1, formula = ~ . + log(Trt), newdata = new_data,
               sample_prior = FALSE, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  expect_error(update(fit1, formula. = ~ . + wrong_var),
               "New variables found: wrong_var")
  up <- update(fit2, algorithm = "fullrank", testmode = TRUE)
  expect_equal(up$algorithm, "fullrank")
  up <- update(fit2, formula. = bf(. ~ ., a + b ~ 1, nl = TRUE), 
               testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  up <- update(fit2, formula. = count ~ a + b, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  
  # VarCorr
  vc <- VarCorr(fit1)
  expect_equal(names(vc), c("visit"))
  Names <- c("Intercept", "Trt")
  expect_equivalent(dimnames(vc$visit$cov$mean), 
                    list(Names, Names))
  expect_output(print(vc), "visit")
  data_vc <- as.data.frame(vc)
  expect_equal(dim(data_vc), c(2, 7))
  expect_equal(names(data_vc), c("Estimate", "Group", "Name", 
                                 "Std.Dev", "Cor", "Cov", "Cov"))
  vc <- VarCorr(fit2)
  expect_equal(names(vc), c("patient"))
  data_vc <- as.data.frame(vc)
  expect_equal(dim(data_vc), c(2, 7))
  
  # vcov
  expect_equal(dim(vcov(fit1)), c(8, 8))
  expect_equal(dim(vcov(fit1, cor = TRUE)), c(8, 8))
  
  # WAIC
  waic1 <- SW(WAIC(fit1))
  expect_true(is.numeric(waic1[["waic"]]))
  expect_true(is.numeric(waic1[["se_waic"]]))
  expect_equal(waic1, SW(waic(fit1)))
  
  waic_compare <- SW(WAIC(fit1, fit1))
  expect_equal(length(waic_compare), 2)
  expect_equal(dim(attr(waic_compare, "compare")), c(1,2))
  waic2 <- SW(WAIC(fit2))
  expect_true(is.numeric(waic2[["waic"]]))
  waic_pointwise <- SW(WAIC(fit2, pointwise = TRUE))
  expect_equal(waic2, waic_pointwise)
  expect_warning(WAIC(fit1, fit2), "Model comparisons are most likely invalid")
  waic4 <- SW(WAIC(fit4))
  expect_true(is.numeric(waic4[["waic"]]))
  
  # test diagnostic convenience functions
  expect_true(is(log_posterior(fit1), "data.frame"))
  expect_true(is(nuts_params(fit1), "data.frame"))
  expect_true(is(rhat(fit1), "numeric"))
  expect_true(is(neff_ratio(fit1), "numeric"))
})