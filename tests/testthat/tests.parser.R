test_that("Test that extract.effects finds all variables", {
  expect_equal(extract.effects(t2_brand_recall ~ psi_expsi + psi_api_probsolv + 
                                 psi_api_ident + psi_api_intere + psi_api_groupint)$all, 
               t2_brand_recall ~ psi_expsi + psi_api_probsolv + psi_api_ident + psi_api_intere + psi_api_groupint)
})

test_that("Test that extract.effects handles addition arguments correctly", {
  expect_equal(extract.effects(y | se(I(a+2)) ~ x, family = "gaussian")$se, ~I(a+2))
  expect_equal(extract.effects(y | se(I(a+2)) ~ x, family = "gaussian")$all, y~ x + I(a+2))
  expect_equal(extract.effects(y | weights(I(1/n)) ~ x, family = "gaussian")$weights, ~I(1/n))
  expect_equal(extract.effects(y | se(I(a+2)) | cens(log(b)) ~ x, family = "gaussian")$cens, ~log(b))
  expect_equal(extract.effects(y | trials(10) ~ x, family = "binomial")$trials, 10)
  expect_equal(extract.effects(y | cat(cate) ~ x, family = "cumulative")$cat, ~cate)
  expect_equal(extract.effects(y | cens(I(cens^2)) ~ z + (x|patient), family = "weibull")$all, 
               y ~ z + x + patient + I(cens^2))
})

test_that("Test that brm.update.formula returns correct formulas", {
  expect_equal(brm.update.formula(y~x, addition = list(se = ~I(sei+2))), y | se(I(sei+2)) ~ x)
  expect_equal(brm.update.formula(y~x, addition = list(se = ~sei, cens = ~censored)), 
               y | se(sei) | cens(censored) ~ x)
})

test_that("Test that is.formula is TRUE for formulas and otherwise FALSE", {
  expect_equal(is.formula(y~1), TRUE)
  expect_equal(is.formula("a"), FALSE)
  expect_equal(is.formula(list(y~1, ~1)), TRUE)
  expect_equal(is.formula(list(y~1,1)), TRUE)
  expect_equal(is.formula(list("a",1)), FALSE)
  expect_equal(is.formula(list(y~1, ~1), or = FALSE), TRUE)
  expect_equal(is.formula(list(y~1,1), or = FALSE), FALSE)
})