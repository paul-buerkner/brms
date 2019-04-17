context("Tests for plotting functions")

test_that("plotting functions don't throw unexpected errors", {
  ggplot2::theme_set(theme_black())
  fit1 <- brms:::rename_pars(brms:::brmsfit_example1)
  fit3 <- brms:::rename_pars(brms:::brmsfit_example3)
  
  # plot.brmsfit
  expect_silent(p <- plot(fit1, plot = FALSE))
  expect_silent(p <- plot(fit1, pars = "^b", plot = FALSE))
  expect_silent(p <- plot(fit1, pars = "^sd", plot = FALSE))
  expect_error(plot(fit1, pars = "123"),  "No valid parameters selected")
  
  # stanplot.brmsfit
  expect_silent(p <- stanplot(fit1))
  expect_silent(p <- stanplot(fit1, pars = "^b"))
  expect_silent(p <- suppressMessages(
    stanplot(fit1, type = "trace", pars = "^b_")
  ))
  expect_silent(p <- stanplot(fit1, type = "hist", pars = "^sd_"))
  expect_silent(p <- stanplot(fit1, type = "dens"))
  expect_silent(p <- stanplot(fit1, type = "scatter",
                              pars = parnames(fit1)[2:3], 
                              exact_match = TRUE))
  expect_silent(p <- stanplot(fit1, type = "rhat", pars = "^b_"))
  expect_silent(p <- stanplot(fit1, type = "neff"))
  expect_silent(p <- stanplot(fit1, type = "acf"))
  expect_silent(p <- stanplot(fit1, type = "nuts_divergence"))
  expect_error(stanplot(fit1, type = "density"), "Invalid plot type")
  expect_error(stanplot(fit1, type = "hex"), 
               "Exactly 2 parameters must be selected")
  
  # pairs.brmsfit
  expect_s3_class(SW(pairs(fit1, pars = parnames(fit1)[1:3])), "bayesplot_grid")
  
  # marginal_effects: manual checks of plotting method
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
  expect_true(is(marg_plot[[1]], "ggplot"))
  # test with 1 categorical predictor
  attr(marg_results[[1]], "effects") <- "P2"
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
  # test with 1 numeric and 1 categorical predictor
  attr(marg_results[[1]], "effects") <- c("P1", "P2")
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
  # test ordinal raster plot
  attr(marg_results[[1]], "effects") <- c("P1", "cats__")
  attr(marg_results[[1]], "ordinal") <- TRUE
  marg_plot <- plot(marg_results, plot = FALSE)
  expect_true(is(marg_plot[[1]], "ggplot"))
})
