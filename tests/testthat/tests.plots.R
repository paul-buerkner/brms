context("Tests for plotting functions")

test_that("plotting functions don't throw unexpected errors", {
  theme_set(theme_black())
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  
  # plot.brmsfit
  expect_silent(p <- plot(fit, plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^b", plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^sd", plot = FALSE))
  expect_error(plot(fit, pars = "123"),  "No valid parameters selected")
  
  # stanplot.brmsfit
  expect_silent(p <- stanplot(fit))
  expect_silent(p <- stanplot(fit, pars = "^b"))
  expect_silent(p <- suppressMessages(
    stanplot(fit, type = "trace", pars = "^b_")
  ))
  expect_silent(p <- stanplot(fit, type = "hist", pars = "^sd_"))
  expect_silent(p <- stanplot(fit, type = "dens"))
  expect_silent(p <- stanplot(fit, type = "scatter",
                              pars = parnames(fit)[2:3], 
                              exact_match = TRUE))
  expect_silent(p <- stanplot(fit, type = "rhat", pars = "^b_"))
  expect_silent(p <- stanplot(fit, type = "neff"))
  expect_silent(p <- stanplot(fit, type = "acf"))
  expect_silent(p <- stanplot(fit, type = "nuts_divergence"))
  expect_error(stanplot(fit, type = "density"), "Invalid plot type")
  expect_error(stanplot(fit, type = "hex"), 
               "Exactly 2 parameters must be selected")
  
  # pairs.brmsfit
  expect_s3_class(SW(pairs(fit, pars = parnames(fit)[1:3])), "bayesplot_grid")
  
  # marginal_effects: manual checks of plotting method
  N <- 90
  marg_results <- data.frame(
    P1 = rpois(N, 20), 
    P2 = factor(rep(1:3, each = N / 3)),
    estimate__ = rnorm(N, sd = 5), 
    se__ = rt(N, df = 10), 
    cond__ = rep(1:2, each = N / 2),
    cats__ = factor(rep(1:3, each = N / 3))
  )
  marg_results[["lower__"]] <- marg_results$estimate__ - 2
  marg_results[["upper__"]] <- marg_results$estimate__ + 2
  marg_results <- list(marg_results[order(marg_results$P1), ])
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
