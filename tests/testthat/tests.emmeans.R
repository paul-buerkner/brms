context("Tests for emmeans support")

skip_on_cran()

require(emmeans)

SW <- suppressWarnings
fit1 <- rename_pars(brms:::brmsfit_example1)
fit2 <- rename_pars(brms:::brmsfit_example2)
fit4 <- rename_pars(brms:::brmsfit_example4)
fit6 <- rename_pars(brms:::brmsfit_example6)

test_that("emmeans returns expected output structure", {
  em <- summary(emmeans(fit1, "Age", by = "Trt"))
  expect_equal(nrow(em), 2)

  em <- summary(emmeans(fit1, "Trt", dpar = "sigma"))
  expect_equal(nrow(em), 2)

  em <- summary(emmeans(fit1, "Age", by = "Exp"))
  expect_equal(nrow(em), 5)

  em <- summary(emmeans(fit1, "Exp"))
  expect_equal(nrow(em), 5)

  em <- SW(summary(emmeans(fit2, "Age", nlpar = "a")))
  expect_equal(nrow(em), 1)

  em <- SW(summary(emmeans(fit4, "x1", dpar = "mu")))
  expect_equal(nrow(em), 1)
})

test_that("emmeans supports 'epred' predictions", {
  em <- summary(emmeans(fit2, "Age", epred = TRUE))
  expect_equal(nrow(em), 1)

  em <- summary(emmeans(fit2, "Age", by = "Trt", epred = TRUE))
  expect_equal(nrow(em), 2)

  # test for a multivariate model
  em <- summary(emmeans(fit6, "Age", by = "Trt", epred = TRUE))
  expect_equal(nrow(em), 2)
})

test_that("emmeans supports multilevel terms", {
  em <- summary(emmeans(fit1, "Age", by = "Trt", re_formula = NULL))
  expect_equal(nrow(em), 2)

  em <- SW(summary(emmeans(fit2, "Age", nlpar = "a", re_formula = NULL)))
  expect_equal(nrow(em), 1)
})
