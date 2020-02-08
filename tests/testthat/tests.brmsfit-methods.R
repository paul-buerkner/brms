context("Tests for brmsfit methods")

expect_range <- function(object, lower = -Inf, upper = Inf, ...) {
  testthat::expect_true(all(object >= lower & object <= upper), ...)
}
expect_ggplot <- function(object, ...) {
  testthat::expect_true(is(object, "ggplot"), ...)
}

SM <- suppressMessages
SW <- suppressWarnings

fit1 <- brms:::rename_pars(brms:::brmsfit_example1)
fit2 <- brms:::rename_pars(brms:::brmsfit_example2)
fit3 <- brms:::rename_pars(brms:::brmsfit_example3)
fit4 <- brms:::rename_pars(brms:::brmsfit_example4)
fit5 <- brms:::rename_pars(brms:::brmsfit_example5)
fit6 <- brms:::rename_pars(brms:::brmsfit_example6)

# test S3 methods in alphabetical order
test_that("as.data.frame has reasonable ouputs", {
  ps <- as.data.frame(fit1)
  expect_true(is(ps, "data.frame"))
  expect_equal(dim(ps), c(nsamples(fit1), length(parnames(fit1))))
})

test_that("as.matrix has reasonable ouputs", {
  ps <- as.matrix(fit1)
  expect_true(is(ps, "matrix"))
  expect_equal(dim(ps), c(nsamples(fit1), length(parnames(fit1))))
})

test_that("as.array has reasonable ouputs", {
  ps <- as.array(fit1)
  expect_true(is.array(ps))
  chains <- fit1$fit@sim$chains
  ps_dim <- c(nsamples(fit1) / chains, chains, length(parnames(fit1)))
  expect_equal(dim(ps), ps_dim)
})

test_that("as.mcmc has reasonable ouputs", {
  chains <- fit1$fit@sim$chains
  mc <- as.mcmc(fit1)
  expect_equal(length(mc), chains)
  expect_equal(dim(mc[[1]]), c(nsamples(fit1) / chains, length(parnames(fit1))))
  mc <- as.mcmc(fit1, combine_chains = TRUE)
  expect_equal(dim(mc), c(nsamples(fit1), length(parnames(fit1))))
  # test assumes thin = 1
  expect_equal(dim(as.mcmc(fit1, inc_warmup = TRUE)[[1]]), 
               c(fit1$fit@sim$iter, length(parnames(fit1))))
})

test_that("autocor has reasonable ouputs", {
  expect_true(is.null(SW(autocor(fit1))))
  expect_true(is.null(SW(autocor(fit6, resp = "count"))))
})

test_that("bayes_R2 has reasonable ouputs", {
  fit1 <- add_criterion(fit1, "bayes_R2")
  R2 <- bayes_R2(fit1, summary = FALSE)
  expect_equal(dim(R2), c(nsamples(fit1), 1))
  R2 <- bayes_R2(fit2, newdata = model.frame(fit2)[1:5, ], re_formula = NA)
  expect_equal(dim(R2), c(1, 4))
  R2 <- bayes_R2(fit6)
  expect_equal(dim(R2), c(2, 4))
})

test_that("bayes_factor has reasonable ouputs", {
  # don't test for now as it requires calling Stan's C++ code
})

test_that("bridge_sampler has reasonable ouputs", {
  # don't test for now as it requires calling Stan's C++ code
})

test_that("coef has reasonable ouputs", {
  coef1 <- SM(coef(fit1))
  expect_equal(dim(coef1$visit), c(4, 4, 8))
  coef1 <- SM(coef(fit1, summary = FALSE))
  expect_equal(dim(coef1$visit), c(nsamples(fit1), 4, 8))
  coef2 <- SM(coef(fit2))
  expect_equal(dim(coef2$patient), c(59, 4, 4))
  coef4 <- SM(coef(fit4))
  expect_equal(dim(coef4$subject), c(10, 4, 8))
})

test_that("combine_models has reasonable ouputs", {
  expect_equal(nsamples(combine_models(fit1, fit1)), nsamples(fit1) * 2)
})

test_that("conditional_effects has reasonable ouputs", {
  me <- conditional_effects(fit1, resp = "count")
  expect_equal(nrow(me[[2]]), 100)
  meplot <- plot(me, points = TRUE, rug = TRUE, 
                 ask = FALSE, plot = FALSE)
  expect_ggplot(meplot[[1]])
  
  me <- conditional_effects(fit1, "Trt", select_points = 0.1)
  expect_lt(nrow(attr(me[[1]], "points")), nobs(fit1))
  
  me <- conditional_effects(fit1, "Exp:Age", surface = TRUE, 
                            resolution = 15, too_far = 0.2)
  meplot <- plot(me, plot = FALSE)
  expect_ggplot(meplot[[1]])
  meplot <- plot(me, stype = "raster", plot = FALSE)
  expect_ggplot(meplot[[1]])
  
  me <- conditional_effects(fit1, "Age", spaghetti = TRUE, nsamples = 10)
  expect_equal(nrow(attr(me$Age, "spaghetti")), 1000)
  meplot <- plot(me, plot = FALSE)
  expect_ggplot(meplot[[1]])
  expect_error(
    conditional_effects(fit1, "Age", spaghetti = TRUE, surface = TRUE),
    "Cannot use 'spaghetti' and 'surface' at the same time"
  )
  
  mdata = data.frame(
    Age = c(-0.3, 0, 0.3), 
    count = c(10, 20, 30), 
    Exp = c(1, 3, 5)
  )
  exp_nrow <- nrow(mdata) * 100
  me <- conditional_effects(fit1, effects = "Age", conditions = mdata)
  expect_equal(nrow(me[[1]]), exp_nrow)
  
  mdata$visit <- 1:3
  me <- conditional_effects(fit1, re_formula = NULL, conditions = mdata)
  expect_equal(nrow(me$Age), exp_nrow)
  
  me <- conditional_effects(
    fit1, "Age:Trt", int_conditions = list(Age = rnorm(5))
  )
  expect_equal(nrow(me[[1]]), 200)
  me <- conditional_effects(
    fit1, "Age:Trt", int_conditions = list(Age = quantile)
  )
  expect_equal(nrow(me[[1]]), 200)
  
  expect_error(conditional_effects(fit1, effects = "Trtc"), 
               "All specified effects are invalid for this model")
  expect_warning(conditional_effects(fit1, effects = c("Trtc", "Trt")), 
                 "Some specified effects are invalid for this model")
  expect_error(conditional_effects(fit1, effects = "Trtc:a:b"), 
               "please use the 'conditions' argument")
  
  mdata$visit <- NULL
  mdata$patient <- 1
  expect_equal(nrow(conditional_effects(fit2)[[2]]), 100)
  me <- conditional_effects(fit2, re_formula = NULL, conditions = mdata)
  expect_equal(nrow(me$Age), exp_nrow)
  
  expect_warning(
    me4 <- conditional_effects(fit4),
    "Predictions are treated as continuous variables"
  )
  expect_true(is(me4, "brms_conditional_effects"))
  me4 <- conditional_effects(fit4, "x2", categorical = TRUE)
  expect_true(is(me4, "brms_conditional_effects"))
  
  me5 <- conditional_effects(fit5)
  expect_true(is(me5, "brms_conditional_effects"))
  
  me6 <- conditional_effects(fit6, nsamples = 40)
  expect_true(is(me6, "brms_conditional_effects"))
})

test_that("plot of conditional_effects has reasonable outputs", {
  ggplot2::theme_set(theme_black())
  N <- 90
  marg_results <- data.frame(
    effect1__ = rpois(N, 20), 
    effect2__ = factor(rep(1:3, each = N / 3)),
    estimate__ = rnorm(N, sd = 5), 
    se__ = rt(N, df = 10), 
    cond__ = rep(1:2, each = N / 2),
    cats__ = factor(rep(1:3, each = N / 3))
  )
  marg_results[["lower__"]] <- marg_results$estimate__ - 2
  marg_results[["upper__"]] <- marg_results$estimate__ + 2
  marg_results <- list(marg_results[order(marg_results$effect1__), ])
  class(marg_results) <- "brmsMarginalEffects"
  attr(marg_results[[1]], "response") <- "count"
  # test with 1 numeric predictor
  attr(marg_results[[1]], "effects") <- "P1"
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_ggplot(marg_plot[[1]])
  # test with 1 categorical predictor
  attr(marg_results[[1]], "effects") <- "P2"
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_ggplot(marg_plot[[1]])
  # test with 1 numeric and 1 categorical predictor
  attr(marg_results[[1]], "effects") <- c("P1", "P2")
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_ggplot(marg_plot[[1]])
  # test ordinal raster plot
  attr(marg_results[[1]], "effects") <- c("P1", "cats__")
  attr(marg_results[[1]], "ordinal") <- TRUE
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_ggplot(marg_plot[[1]])
})

test_that("conditional_smooths has reasonable ouputs", {
  ms <- conditional_smooths(fit1)
  expect_equal(nrow(ms[[1]]), 100)
  expect_true(is(ms, "brms_conditional_effects"))
  
  ms <- conditional_smooths(fit1, spaghetti = TRUE, nsamples = 10)
  expect_equal(nrow(attr(ms[[1]], "spaghetti")), 1000)
  
  expect_error(conditional_smooths(fit1, smooths = "s3"),
               "No valid smooth terms found in the model")
  expect_error(conditional_smooths(fit2),
               "No valid smooth terms found in the model")
})

test_that("family has reasonable ouputs", {
  expect_is(family(fit1), "brmsfamily")
  expect_is(family(fit6, resp = "count"), "brmsfamily")
  expect_output(print(family(fit1), links = TRUE), "student.*log.*logm1")
  expect_output(print(family(fit5)), "Mixture.*gaussian.*exponential")
})

test_that("fitted has reasonable outputs", {
  skip_on_cran()
  
  fi <- fitted(fit1)
  expect_equal(dim(fi), c(nobs(fit1), 4))
  expect_equal(colnames(fi), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  
  newdata <- data.frame(
    Age = c(0, -0.2), visit = c(1, 4), Trt = c(0, 1), 
    count = c(20, 13), patient = c(1, 42), Exp = c(2, 4)
  )
  fi <- fitted(fit1, newdata = newdata)
  expect_equal(dim(fi), c(2, 4))
  newdata$visit <- c(1, 6)
  fi <- fitted(fit1, newdata = newdata, 
               allow_new_levels = TRUE)
  expect_equal(dim(fi), c(2, 4))
  
  # fitted values with new_levels
  newdata <- data.frame(
    Age = 0, visit = paste0("a", 1:100), Trt = 0, 
    count = 20, patient = 1, Exp = 2
  )
  fi <- fitted(fit1, newdata = newdata, allow_new_levels = TRUE, 
               sample_new_levels = "old_levels", nsamples = 10)
  expect_equal(dim(fi), c(100, 4))
  fi <- fitted(fit1, newdata = newdata, allow_new_levels = TRUE, 
               sample_new_levels = "gaussian", nsamples = 1)
  expect_equal(dim(fi), c(100, 4))
  
  # fitted values of auxiliary parameters
  newdata <- data.frame(
    Age = 0, visit = c("a", "b"), Trt = 0,
    count = 20, patient = 1, Exp = 2
  )
  fi <- fitted(fit1, dpar = "sigma")
  expect_equal(dim(fi), c(nobs(fit1), 4))
  expect_true(all(fi > 0))
  fi_lin <- fitted(fit1, dpar = "sigma", scale = "linear")
  expect_equal(dim(fi_lin), c(nobs(fit1), 4))
  expect_true(!isTRUE(all.equal(fi, fi_lin)))
  expect_error(fitted(fit1, dpar = "inv"),
               "Invalid argument 'dpar'")
  
  fi <- fitted(fit2)
  expect_equal(dim(fi), c(nobs(fit2), 4))
  fi <- fitted(fit2, newdata = newdata,
               allow_new_levels = TRUE)
  expect_equal(dim(fi), c(2, 4))
  fi <- fitted(fit2, dpar = "shape")
  expect_equal(dim(fi), c(nobs(fit2), 4))
  expect_equal(fi[1, ], fi[2, ])
  fi <- fitted(fit2, nlpar = "a")
  expect_equal(dim(fi), c(nobs(fit2), 4))
  
  fi <- fitted(fit3, newdata = fit3$data[1:10, ])
  expect_equal(dim(fi), c(10, 4))
  
  fi <- fitted(fit4)
  expect_equal(dim(fi), c(nobs(fit4), 4, 4))
  fi <- fitted(fit4, newdata = fit4$data[1, ])
  expect_equal(dim(fi), c(1, 4, 4))
  fi <- fitted(fit4, newdata = fit4$data[1, ], scale = "linear")
  expect_equal(dim(fi), c(1, 4, 3))
  
  fi <- fitted(fit5)
  expect_equal(dim(fi), c(nobs(fit5), 4))
  
  fi <- fitted(fit6)
  expect_equal(dim(fi), c(nobs(fit6), 4, 2))
  expect_equal(dimnames(fi)[[3]], c("volume", "count"))
})

test_that("fixef has reasonable ouputs", {
  fixef1 <- SM(fixef(fit1))
  expect_equal(rownames(fixef1), 
               c("Intercept", "sigma_Intercept", "Trt1", "Age", 
                 "Trt1:Age", "sigma_Trt1", "sAge_1", "moExp")
  )
  fixef1 <- SM(fixef(fit1, pars = c("Age", "sAge_1")))
  expect_equal(rownames(fixef1), c("Age", "sAge_1"))
})

test_that("formula has reasonable ouputs", {
  expect_true(is.brmsformula(formula(fit1))) 
})

test_that("hypothesis has reasonable ouputs", {
  hyp <- hypothesis(fit1, c("Age > Trt1", "Trt1:Age = -1"))
  expect_equal(dim(hyp$hypothesis), c(2, 8))
  expect_output(print(hyp), "(Age)-(Trt1) > 0", fixed = TRUE)
  expect_ggplot(plot(hyp, plot = FALSE)[[1]])
  
  hyp <- hypothesis(fit1, "Intercept = 0", class = "sd", group = "visit")
  expect_true(is.numeric(hyp$hypothesis$Evid.Ratio[1]))
  expect_output(print(hyp), "class sd_visit:", fixed = TRUE)
  expect_ggplot(plot(hyp, ignore_prior = TRUE, plot = FALSE)[[1]])
  
  hyp <- hypothesis(fit1, "0 > r_visit[4,Intercept]", class = "", alpha = 0.01)
  expect_equal(dim(hyp$hypothesis), c(1, 8))
  expect_output(print(hyp, chars = NULL), "r_visit[4,Intercept]", fixed = TRUE)
  expect_output(print(hyp), "99%-CI", fixed = TRUE)
  
  hyp <- hypothesis(
    fit1, c("Intercept = 0", "Intercept + exp(Trt1) = 0"),
    group = "visit", scope = "coef"
  )
  expect_equal(dim(hyp$hypothesis), c(8, 9))
  expect_equal(hyp$hypothesis$Group[1], factor(1, levels = 1:4))
  
  expect_error(hypothesis(fit1, "Intercept > x"), fixed = TRUE,
               "cannot be found in the model: \n'b_x'")
  expect_error(hypothesis(fit1, 1),
               "Argument 'hypothesis' must be a character vector")
  expect_error(hypothesis(fit2, "b_Age = 0", alpha = 2),
               "Argument 'alpha' must be a single value in [0,1]",
               fixed = TRUE)
  expect_error(hypothesis(fit3, "b_Age x 0"),
               "Every hypothesis must be of the form 'left (= OR < OR >) right'",
               fixed = TRUE)
  
  # test hypothesis.default method
  hyp <- hypothesis(as.data.frame(fit3), "bsp_meAgeAgeSD > sigma")
  expect_equal(dim(hyp$hypothesis), c(1, 8))
  hyp <- hypothesis(fit3$fit, "bsp_meAgeAgeSD > sigma")
  expect_equal(dim(hyp$hypothesis), c(1, 8))
})

test_that("launch_shinystan has reasonable ouputs", {
  # requires running shiny which is not reasonable in automated tests
})

test_that("log_lik has reasonable ouputs", {
  expect_equal(dim(log_lik(fit1)), c(nsamples(fit1), nobs(fit1)))
  expect_equal(dim(logLik(fit1)), c(nsamples(fit1), nobs(fit1)))
  expect_equal(dim(log_lik(fit2)), c(nsamples(fit2), nobs(fit2)))
})

test_that("loo has reasonable outputs", {
  skip_on_cran()
  
  loo1 <- SW(LOO(fit1, cores = 1))
  expect_true(is.numeric(loo1$estimates))
  expect_output(print(loo1), "looic")
  
  loo_compare1 <- SW(loo(fit1, fit1, cores = 1))
  expect_equal(names(loo_compare1$loos), c("fit1", "fit1"))
  expect_equal(dim(loo_compare1$ic_diffs__), c(1, 2))
  expect_output(print(loo_compare1), "'fit1':")
  expect_is(loo_compare1$diffs, "compare.loo")
  
  loo2 <- SW(loo(fit2, cores = 1))
  expect_true(is.numeric(loo2$estimates))
  
  loo3 <- SW(loo(fit3, cores = 1))
  expect_true(is.numeric(loo3$estimates))
  loo3 <- SW(loo(fit3, pointwise = TRUE, cores = 1))
  expect_true(is.numeric(loo3$estimates))
  
  loo4 <- SW(loo(fit4, cores = 1))
  expect_true(is.numeric(loo4$estimates))
  
  # fails because of too small effective sample size
  # loo5 <- SW(loo(fit5, cores = 1))
  # expect_true(is.numeric(loo5$estimates))
  
  loo6_1 <- SW(loo(fit6, cores = 1))
  expect_true(is.numeric(loo6_1$estimates))
  loo6_2 <- SW(loo(fit6, cores = 1, newdata = fit6$data))
  expect_true(is.numeric(loo6_2$estimates))
  loo_compare <- loo_compare(loo6_1, loo6_2)
  expect_range(loo_compare[2, 1], -1, 1)
})

test_that("loo_subsample has reasonable outputs", {
  skip_on_cran()
  
  loo2 <- SW(loo_subsample(fit2, observations = 50))
  expect_true(is.numeric(loo2$estimates))
  expect_equal(nrow(loo2$pointwise), 50)
  expect_output(print(loo2), "looic")
})

test_that("loo_R2 has reasonable outputs", {
  skip_on_cran()
  
  R2 <- SW(loo_R2(fit1))
  expect_equal(length(R2), 1)
  
  # fails on travis for some strange reason
  # R2 <- SW(loo_R2(fit6))
  # expect_equal(length(R2), 2)
})

test_that("loo_linpred has reasonable outputs", {
  skip_on_cran()
  
  llp <- SW(loo_linpred(fit1))
  expect_equal(length(llp), nobs(fit1))
  expect_error(loo_linpred(fit4), "Method 'loo_linpred'")
  llp <- SW(loo_linpred(fit2, scale = "response", type = "var"))
  expect_equal(length(llp), nobs(fit2))
})

test_that("loo_predict has reasonable outputs", {
  skip_on_cran()
  
  llp <- SW(loo_predict(fit1))
  expect_equal(length(llp), nobs(fit1))
  
  newdata <- data.frame(
    Age = 0, visit = c("a", "b"), Trt = 0,
    count = 20, patient = 1, Exp = 2
  )
  llp <- SW(loo_predict(
    fit1, newdata = newdata,
    type = "quantile", probs = c(0.25, 0.75),
    allow_new_levels = TRUE
  ))
  expect_equal(dim(llp), c(2, nrow(newdata)))
  llp <- SW(loo_predict(fit4))
  expect_equal(length(llp), nobs(fit4))
})

test_that("loo_predictive_interval has reasonable outputs", {
  skip_on_cran()
  
  llp <- SW(loo_predictive_interval(fit3))
  expect_equal(dim(llp), c(nobs(fit3), 2))
})

test_that("loo_model_weights has reasonable outputs", {
  skip_on_cran()
  
  llw <- SW(loo_model_weights(fit1, fit1))
  expect_is(llw[1:2], "numeric")
  expect_equal(names(llw), c("fit1", "fit1"))
})

test_that("model.frame has reasonable ouputs", {
  expect_equal(model.frame(fit1), fit1$data)
})

test_that("model_weights has reasonable ouputs", {
  mw <- model_weights(fit1, fit1, weights = "waic")
  expect_equal(mw, setNames(c(0.5, 0.5), c("fit1", "fit1")))
})

test_that("ngrps has reasonable ouputs", {
  expect_equal(ngrps(fit1), list(visit = 4))
  expect_equal(ngrps(fit2), list(patient = 59))
})

test_that("nobs has reasonable ouputs", {
  expect_equal(nobs(fit1), nrow(epilepsy))
})

test_that("nsamples has reasonable ouputs", {
  expect_equal(nsamples(fit1), 50)
  expect_equal(nsamples(fit1, subset = 10:1), 10)
  expect_equal(nsamples(fit1, incl_warmup = TRUE), 200)
})

test_that("pairs has reasonable outputs", {
  expect_s3_class(SW(pairs(fit1, pars = parnames(fit1)[1:3])), 
                  "bayesplot_grid")
})

test_that("parnames has reasonable ouputs", {
  expect_true(all(
    c("b_Intercept", "bsp_moExp", "ar[1]", "cor_visit__Intercept__Trt1", 
      "nu", "simo_moExp1[2]", "r_visit[4,Trt1]", "s_sAge_1[8]", 
      "prior_sd_visit", "prior_cor_visit", "lp__") %in%
      parnames(fit1)  
  ))
  expect_true(all(
    c("b_a_Intercept", "b_b_Age", "sd_patient__b_Intercept",
      "cor_patient__a_Intercept__b_Intercept", 
      "r_patient__a[1,Intercept]", "r_patient__b[4,Intercept]",
      "prior_b_a") %in%
      parnames(fit2)  
  ))
  expect_true(all(
    c("lscale_volume_gpAgeTrt0", "lscale_volume_gpAgeTrt1") %in% 
      parnames(fit6)
  ))
})

test_that("plot has reasonable outputs", {
  expect_silent(p <- plot(fit1, plot = FALSE))
  expect_silent(p <- plot(fit1, pars = "^b", plot = FALSE))
  expect_silent(p <- plot(fit1, pars = "^sd", plot = FALSE))
  expect_error(plot(fit1, pars = "123"),  "No valid parameters selected")
})

test_that("post_prob has reasonable ouputs", {
  # only test error messages for now
  expect_error(post_prob(fit1, fit2, model_names = "test1"),
               "Number of model names is not equal to the number of models")
})

test_that("posterior_average has reasonable outputs", {
  pnames <- c("b_Age", "nu")
  ps <- posterior_average(fit1, fit1, pars = pnames, weights = c(0.3, 0.7))
  expect_equal(dim(ps), c(nsamples(fit1), 2))
  expect_equal(names(ps), pnames)
  
  weights <- rexp(3)
  ps <- brms:::SW(posterior_average(
    fit1, fit2, fit3, pars = "nu", weights = rexp(3), 
    missing = 1, nsamples = 10
  ))
  expect_equal(dim(ps), c(10, 1))
  expect_equal(names(ps), "nu")
})

test_that("posterior_samples has reasonable outputs", {
  ps <- posterior_samples(fit1)
  expect_equal(dim(ps), c(nsamples(fit1), length(parnames(fit1))))
  expect_equal(names(ps), parnames(fit1))
  expect_equal(names(posterior_samples(fit1, pars = "^b_")),
               c("b_Intercept", "b_sigma_Intercept", "b_Trt1", 
                 "b_Age", "b_Trt1:Age", "b_sigma_Trt1"))
  
  # test default method
  ps <- posterior_samples(fit1$fit, "^b_Intercept$")
  expect_equal(dim(ps), c(nsamples(fit1), 1))
})

test_that("posterior_summary has reasonable outputs", {
  ps <- posterior_summary(fit1, "^b_")
  expect_equal(dim(ps), c(6, 4))
})

test_that("posterior_interval has reasonable outputs", {
  expect_equal(dim(posterior_interval(fit1)), 
               c(length(parnames(fit1)), 2))
})

test_that("posterior_predict has reasonable outputs", {
  expect_equal(dim(posterior_predict(fit1)), 
               c(nsamples(fit1), nobs(fit1)))
})

test_that("posterior_linpred has reasonable outputs", {
  expect_equal(dim(posterior_linpred(fit1)), 
               c(nsamples(fit1), nobs(fit1)))
})

test_that("pp_average has reasonable outputs", {
  ppa <- pp_average(fit1, fit1, weights = "waic")
  expect_equal(dim(ppa), c(nobs(fit1), 4))
  ppa <- pp_average(fit1, fit1, weights = c(1, 4))
  expect_equal(attr(ppa, "weights"), c(fit1 = 0.2, fit1 = 0.8))
  ns <- c(fit1 = nsamples(fit1) / 5, fit1 = 4 * nsamples(fit1) / 5)
  expect_equal(attr(ppa, "nsamples"), ns)
})

test_that("pp_check has reasonable outputs", {
  expect_ggplot(pp_check(fit1))
  expect_ggplot(pp_check(fit1, newdata = fit1$data[1:100, ]))
  expect_ggplot(pp_check(fit1, "stat", nsamples = 5))
  expect_ggplot(pp_check(fit1, "error_binned"))
  pp <- pp_check(fit1, "ribbon_grouped", group = "visit", x = "Age")
  expect_ggplot(pp)
  pp <- pp_check(fit1, type = "violin_grouped", 
                 group = "visit", newdata = fit1$data[1:100, ])
  expect_ggplot(pp)
  
  pp <- SW(pp_check(fit1, type = "loo_pit", cores = 1))
  expect_ggplot(pp)
  
  expect_ggplot(pp_check(fit3))
  expect_ggplot(pp_check(fit2, "ribbon", x = "Age"))
  expect_error(pp_check(fit2, "ribbon", x = "x"),
               "Variable 'x' could not be found in the data")
  expect_error(pp_check(fit1, "wrong_type"))
  expect_error(pp_check(fit2, "violin_grouped"), "group")
  expect_error(pp_check(fit1, "stat_grouped", group = "g"),
               "Variable 'g' could not be found in the data")
  expect_ggplot(pp_check(fit4))
  expect_ggplot(pp_check(fit5))
  expect_error(pp_check(fit4, "error_binned"),
               "Type 'error_binned' is not available")
})

test_that("pp_expect has reasonable outputs", {
  expect_equal(dim(pp_expect(fit1)), c(nsamples(fit1), nobs(fit1)))
})

test_that("pp_mixture has reasonable outputs", {
  expect_equal(dim(pp_mixture(fit5)), c(nobs(fit5), 4, 2))
  expect_error(pp_mixture(fit1), 
               "Method 'pp_mixture' can only be applied to mixture models"
  )
})

test_that("predict has reasonable outputs", {
  pred <- predict(fit1)
  expect_equal(dim(pred), c(nobs(fit1), 4))
  expect_equal(colnames(pred), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  pred <- predict(fit1, nsamples = 10, probs = c(0.2, 0.5, 0.8))
  expect_equal(dim(pred), c(nobs(fit1), 5))
  
  newdata <- data.frame(
    Age = c(0, -0.2), visit = c(1, 4), Trt = c(1, 0), 
    count = c(2, 10), patient = c(1, 42), Exp = c(1, 2)
  )
  pred <- predict(fit1, newdata = newdata)
  expect_equal(dim(pred), c(2, 4))
  
  newdata$visit <- c(1, 6)
  pred <- predict(fit1, newdata = newdata, allow_new_levels = TRUE)
  expect_equal(dim(pred), c(2, 4))
  
  # predict NA responses in ARMA models
  df <- fit1$data[1:10, ]
  df$count[8:10] <- NA
  pred <- predict(fit1, newdata = df, nsamples = 1)
  expect_true(!anyNA(pred[, "Estimate"]))
  
  pred <- predict(fit2)
  expect_equal(dim(pred), c(nobs(fit2), 4))
  
  pred <- predict(fit2, newdata = newdata, allow_new_levels = TRUE)
  expect_equal(dim(pred), c(2, 4))
  
  # check if grouping factors with a single level are accepted
  newdata$patient <- factor(2)
  pred <- predict(fit2, newdata = newdata)
  expect_equal(dim(pred), c(2, 4))
  
  pred <- predict(fit4)
  expect_equal(dim(pred), c(nobs(fit4), 4))
  expect_equal(colnames(pred), paste0("P(Y = ", 1:4, ")"))
  pred <- predict(fit4, newdata = fit4$data[1, ])
  expect_equal(dim(pred), c(1, 4))
  
  pred <- predict(fit5)
  expect_equal(dim(pred), c(nobs(fit5), 4))
  newdata <- fit5$data[1:10, ]
  newdata$patient <- "a"
  pred <- predict(fit5, newdata, allow_new_levels = TRUE,
                  sample_new_levels = "old_levels")
  expect_equal(dim(pred), c(10, 4))
  pred <- predict(fit5, newdata, allow_new_levels = TRUE,
                  sample_new_levels = "gaussian")
  expect_equal(dim(pred), c(10, 4))
})

test_that("predictive_error has reasonable outputs", {
  expect_equal(dim(predictive_error(fit1)), 
               c(nsamples(fit1), nobs(fit1)))
})

test_that("print has reasonable outputs", {
  expect_output(SW(print(fit1)), "Group-Level Effects:")
})

test_that("prior_samples has reasonable outputs", {
  prs1 <- prior_samples(fit1)
  prior_names <- c(
    "Intercept", "b", "bsp", paste0("simo_moExp1[", 1:4, "]"), 
    "bs", "sds_sAge_1", "b_sigma", "Intercept_sigma", "nu", 
    "sd_visit", "cor_visit"
  )
  expect_equal(colnames(prs1), prior_names)
  
  prs2 <- prior_samples(fit1, pars = "b_Trt1")
  expect_equal(dimnames(prs2), list(as.character(1:nsamples(fit1)), "b_Trt1"))
  expect_equal(sort(prs1$b), sort(prs2$b_Trt))
  
  # test default method
  prs <- prior_samples(fit1$fit, pars = "^sd_visit")
  expect_equal(names(prs), "prior_sd_visit")
})

test_that("prior_summary has reasonable outputs", {
  expect_true(is(prior_summary(fit1), "brmsprior"))
})

test_that("ranef has reasonable outputs", {
  ranef1 <- SM(ranef(fit1))
  expect_equal(dim(ranef1$visit), c(4, 4, 2))
  
  ranef1 <- SM(ranef(fit1, pars = "Trt1"))
  expect_equal(dimnames(ranef1$visit)[[3]], "Trt1")
  
  ranef1 <- SM(ranef(fit1, groups = "a"))
  expect_equal(length(ranef1), 0L)
  
  ranef2 <- SM(ranef(fit2, summary = FALSE))
  expect_equal(dim(ranef2$patient), c(nsamples(fit2), 59, 2))
})

test_that("residuals has reasonable outputs", {
  res1 <- SW(residuals(fit1, type = "pearson", probs = c(0.65)))
  expect_equal(dim(res1), c(nobs(fit1), 3))
  newdata <- cbind(epilepsy[1:10, ], Exp = rep(1:5, 2))
  res2 <- residuals(fit1, newdata = newdata)
  expect_equal(dim(res2), c(10, 4))
  newdata$visit <- rep(1:5, 2)
  
  res3 <- residuals(fit1, newdata = newdata, allow_new_levels = TRUE)
  expect_equal(dim(res3), c(10, 4))
  
  res4 <- residuals(fit2)
  expect_equal(dim(res4), c(nobs(fit2), 4))
  
  expect_error(residuals(fit4), "Predictive errors are not defined")
  
  res6 <- residuals(fit6)
  expect_equal(dim(res6), c(nobs(fit6), 4, 2))
  expect_equal(dimnames(res6)[[3]], c("volume", "count"))
})

test_that("stancode has reasonable outputs", {
  expect_true(is.character(stancode(fit1)))
  expect_output(print(stancode(fit1)), "generated quantities")
})

test_that("standata has reasonable outputs", {
  expect_equal(sort(names(standata(fit1))),
    sort(c("N", "Y",  "Kar", "Kma", "J_lag", "K", "X", "Ksp", "Imo", 
           "Xmo_1", "Jmo", "con_simo_1", "Z_1_1", "Z_1_2", "nb_1", 
           "knots_1", "Zs_1_1", "Ks", "Xs", "offsets", "K_sigma", 
           "X_sigma", "J_1", "N_1", "M_1", "NC_1", "prior_only"))
  )
  expect_equal(sort(names(standata(fit2))),
    sort(c("N", "Y", "weights", "C_1", "K_a", "X_a", "Z_1_a_1",
           "K_b", "X_b", "Z_1_b_2", "J_1", "N_1", "M_1",
           "NC_1", "prior_only"))
  )
})

test_that("mcmc_plot has reasonable outputs", {
  expect_ggplot(mcmc_plot(fit1))
  expect_ggplot(mcmc_plot(fit1, pars = "^b"))
  expect_ggplot(SM(mcmc_plot(fit1, type = "trace", pars = "^b_")))
  expect_ggplot(mcmc_plot(fit1, type = "hist", pars = "^sd_"))
  expect_ggplot(mcmc_plot(fit1, type = "dens"))
  expect_ggplot(mcmc_plot(fit1, type = "scatter",
                          pars = parnames(fit1)[2:3], 
                          fixed = TRUE))
  expect_ggplot(SW(mcmc_plot(fit1, type = "rhat", pars = "^b_")))
  expect_ggplot(SW(mcmc_plot(fit1, type = "neff")))
  expect_ggplot(mcmc_plot(fit1, type = "acf"))
  expect_silent(p <- mcmc_plot(fit1, type = "nuts_divergence"))
  expect_error(mcmc_plot(fit1, type = "density"), "Invalid plot type")
  expect_error(mcmc_plot(fit1, type = "hex"), 
               "Exactly 2 parameters must be selected")
})

test_that("summary has reasonable outputs", {
  summary1 <- SW(summary(fit1, priors = TRUE))
  expect_true(is.numeric(summary1$fixed))
  expect_equal(rownames(summary1$fixed), 
               c("Intercept", "sigma_Intercept", "Trt1", "Age", 
                 "Trt1:Age", "sigma_Trt1", "sAge_1", "moExp"))
  expect_equal(colnames(summary1$fixed), 
               c("Estimate", "Est.Error", "l-95% CI", 
                 "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS"))
  expect_equal(rownames(summary1$random$visit), 
               c("sd(Intercept)", "sd(Trt1)", "cor(Intercept,Trt1)"))
  expect_output(print(summary1), "Population-Level Effects:")
  expect_output(print(summary1), "Priors:")
  
  summary5 <- SW(summary(fit5))
  expect_output(print(summary5), "sigma1")
  expect_output(print(summary5), "theta1")
  
  summary6 <- SW(summary(fit6))
  expect_output(print(summary6), "sdgp")
})

test_that("update has reasonable outputs", {
  # Do not actually refit the model as is causes CRAN checks to fail.
  # Some tests are commented out as they fail when updating Stan code
  # of internal example models because of Stan code mismatches. Refitting
  # these example models is slow especially when done repeatedly and
  # leads the git repo to blow up eventually due the size of the models.
  up <- update(fit1, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  
  new_data <- data.frame(
    Age = rnorm(18), visit = rep(c(3, 2, 4), 6),
    Trt = rep(0:1, 9), count = rep(c(5, 17, 28), 6),
    patient = 1, Exp = 4
  )
  up <- update(fit1, newdata = new_data, save_ranef = FALSE, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  expect_equal(up$data.name, "new_data")
  # expect_equal(attr(up$ranef, "levels")$visit, c("2", "3", "4"))
  # expect_true("r_1_1" %in% up$exclude)
  expect_error(update(fit1, data = new_data), "use argument 'newdata'")
  
  up <- update(fit1, formula = ~ . + I(exp(Age)), testmode = TRUE,
               prior = set_prior("normal(0,10)"))
  expect_true(is(up, "brmsfit"))
  up <- update(fit1, ~ . - Age + factor(Age),  testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  
  up <- update(fit1, formula = ~ . + I(exp(Age)), newdata = new_data,
               sample_prior = FALSE, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  expect_error(update(fit1, formula. = ~ . + wrong_var),
               "New variables found: 'wrong_var'")
  
  up <- update(fit1, save_ranef = FALSE, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  # expect_true("r_1_1" %in% up$exclude)
  up <- update(fit3, save_mevars = FALSE, testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  # expect_true("Xme_1" %in% up$exclude)
  
  up <- update(fit2, algorithm = "fullrank", testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  # expect_equal(up$algorithm, "fullrank")
  up <- update(fit2, formula. = bf(. ~ ., a + b ~ 1, nl = TRUE), 
               testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  up <- update(fit2, formula. = bf(count ~ a + b, nl = TRUE), testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  up <- update(fit3, family = acat(), testmode = TRUE)
  expect_true(is(up, "brmsfit"))
  up <- update(fit3, bf(~., family = acat()), testmode = TRUE)
  expect_true(is(up, "brmsfit"))
})

test_that("VarCorr has reasonable outputs", {
  vc <- VarCorr(fit1)
  expect_equal(names(vc), c("visit"))
  Names <- c("Intercept", "Trt1")
  expect_equal(dimnames(vc$visit$cov)[c(1, 3)], list(Names, Names))
  vc <- VarCorr(fit2)
  expect_equal(names(vc), c("patient"))
  expect_equal(dim(vc$patient$cor), c(2, 4, 2))
  vc <- VarCorr(fit2, summary = FALSE)
  expect_equal(dim(vc$patient$cor), c(nsamples(fit2), 2, 2))
  expect_equal(dim(VarCorr(fit6)$residual__$sd), c(1, 4))
  vc <- VarCorr(fit5)
  expect_equal(dim(vc$patient$sd), c(2, 4))
})

test_that("vcov has reasonable outputs", {
  expect_equal(dim(vcov(fit1)), c(8, 8))
  expect_equal(dim(vcov(fit1, cor = TRUE)), c(8, 8))
})

test_that("waic has reasonable outputs", {
  waic1 <- SW(WAIC(fit1))
  expect_true(is.numeric(waic1$estimates))
  expect_equal(waic1, SW(waic(fit1)))
  
  fit1 <- SW(add_criterion(fit1, "waic"))
  expect_equal(waic(fit1), fit1$criteria$waic)
  
  waic_compare <- SW(waic(fit1, fit1))
  expect_equal(length(waic_compare$loos), 2)
  expect_equal(dim(waic_compare$ic_diffs__), c(1, 2))
  waic2 <- SW(waic(fit2))
  expect_true(is.numeric(waic2$estimates))
  waic_pointwise <- SW(waic(fit2, pointwise = TRUE))
  expect_equal(waic2, waic_pointwise)
  expect_warning(compare_ic(waic1, waic2), 
                 "Model comparisons are likely invalid")
  waic4 <- SW(waic(fit4))
  expect_true(is.numeric(waic4$estimates))
})

test_that("diagnostic convenience functions have reasonable outputs", {
  expect_true(is(log_posterior(fit1), "data.frame"))
  expect_true(is(nuts_params(fit1), "data.frame"))
  expect_true(is(rhat(fit1), "numeric"))
  expect_true(is(neff_ratio(fit1), "numeric"))
})

test_that("contrasts of grouping factors are not stored #214", {
  expect_true(is.null(attr(fit1$data$patient, "contrasts")))
})
