test_that("plotting functions don't throw unexpected errors", {
  skip_if_not(exists("brmsfit_example", asNamespace("brms")))
  fit <- rename_pars(brmsfit_example)
  
  # plot.brmsfit
  expect_silent(p <- plot(fit, do_plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^b", do_plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^sd", do_plot = FALSE))
  expect_error(plot(fit, pars = "123"),  "No valid parameters selected")
  
  # stanplot.brmsfit
  expect_silent(p <- stanplot(fit, quiet = TRUE))
  expect_silent(p <- stanplot(fit, pars = "^b", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "trace", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "hist", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "dens", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "scat", quiet = TRUE,
                              pars = parnames(fit)[2:3], 
                              exact_match = TRUE))
  #expect_silent(p <- stanplot(fit, type = "diag", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "rhat", pars = "^b_",
                              quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "ess", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "mcse", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "ac", quiet = TRUE))
  # warning occurs somewhere in rstan
  expect_silent(SW(stanplot(fit, type = "par", 
                            pars = "^b_Intercept$")))
  expect_error(p <- stanplot(fit, type = "par", pars = "^b_"),
               "'stan_par' expects a single parameter name")
  expect_error(stanplot(fit, type = "density"), "Invalid plot type")
  
  # pairs.brmsfit
  expect_identical(SW(pairs(fit, pars = parnames(fit)[1:3])), NULL)
  
  # marginal_effects
  mdata = data.frame(Age = c(-0.3, 0, 0.3), count = c(10, 20, 30), 
                     visit = 1:3, patient = 1, Trt = 0, Exp = c(1,3,5))
  me <- marginal_effects(fit, conditions = mdata)
  marg_plot <- plot(me, points = TRUE, rug = TRUE, do_plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
  # some manual checks
  N <- 90
  marg_results <- data.frame(P1 = rpois(N, 20), 
                             P2 = factor(rep(1:3, each = N / 3)),
                             Estimate = rnorm(N, sd = 5), 
                             Est.Error = rt(N, df = 10), 
                             MargRow = rep(1:2, each = N / 2))
  marg_results[["lowerCI"]] <- marg_results$Estimate - 2
  marg_results[["upperCI"]] <- marg_results$Estimate + 2
  marg_results <- list(marg_results[order(marg_results$P1), ])
  class(marg_results) <- "brmsMarginalEffects"
  attr(marg_results[[1]], "response") <- "count"
  # test with 1 numeric predictor
  attr(marg_results[[1]], "effects") <- "P1"
  marg_plot <- plot(marg_results, do_plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
  # test with 1 categorical predictor
  attr(marg_results[[1]], "effects") <- "P2"
  marg_plot <- plot(marg_results, do_plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
  # test with 1 numeric and 1 categorical predictor
  attr(marg_results[[1]], "effects") <- c("P1", "P2")
  marg_plot <- plot(marg_results, do_plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
})
