context("Tests for emmeans support")

require(emmeans)
fit1 <- rename_pars(brms:::brmsfit_example1)

test_that("emmeans returns expected output structure", {
  em <- summary(emmeans(fit1, "Age", by = "Trt"))
  expect_equal(nrow(em), 2)
  
  em <- summary(emmeans(fit1, "Trt", dpar = "sigma"))
  expect_equal(nrow(em), 2)
})
