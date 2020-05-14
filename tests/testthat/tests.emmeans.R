context("Tests for emmeans support")

require(emmeans)

SW <- suppressWarnings
fit1 <- rename_pars(brms:::brmsfit_example1)
fit2 <- rename_pars(brms:::brmsfit_example2)
fit4 <- rename_pars(brms:::brmsfit_example4)

test_that("emmeans returns expected output structure", {
  em <- summary(emmeans(fit1, "Age", by = "Trt"))
  expect_equal(nrow(em), 2)
  
  em <- summary(emmeans(fit1, "Trt", dpar = "sigma"))
  expect_equal(nrow(em), 2)
  
  em <- SW(summary(emmeans(fit2, "Age", nlpar = "a")))
  expect_equal(nrow(em), 1)
  
  em <- summary(emmeans(fit4, "x1"))
  expect_equal(nrow(em), 1)
})
