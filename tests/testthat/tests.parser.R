test_that("Test that extract.effects finds all variables in very long formulas", {
  expect_equal(extract.effects(t2_brand_recall ~ psi_expsi + psi_api_probsolv + 
                                 psi_api_ident + psi_api_intere + psi_api_groupint)$all, 
               t2_brand_recall ~ psi_expsi + psi_api_probsolv + psi_api_ident + psi_api_intere + psi_api_groupint)
})

test_that("Test that extract.effects finds all random effects and grouping factors", {
  expect_equal(extract.effects(y ~ a + (1+x|g2) + x + (1|g1) + z)$random, list(~1 + x, ~1))
  expect_equal(extract.effects(y ~ (1+x|g2) + x + (1|g1))$group, list("g2", "g1"))
  expect_equal(extract.effects(y ~ (1+x|g1:g2) + x + (1|g1))$group, list("g1__g2", "g1"))
  expect_error(extract.effects(y ~ (1+x|g1/g2) + x + (1|g1)), 
    paste("Illegal grouping term: g1/g2 \nGrouping terms may contain only variable names",
          "combined by the interaction symbol ':'\n"))
})

test_that("Test that extract.effects accepts || syntax", {
  expect_equal(extract.effects(y ~ a + (1+x||g2) + (1+z|g1))$cor, list(FALSE,TRUE))
  expect_equal(extract.effects(y ~ a + (1+x||g2))$random, list(~1 + x))
  expect_equal(extract.effects(y ~ (1+x||g2) + x + (1||g1))$group, list("g2", "g1"))
  expect_equal(extract.effects(y ~ (1+x||g1:g2))$group, list("g1__g2"))
  expect_error(extract.effects(y ~ (1+x||g1/g2) + x + (1|g1)), 
               paste("Illegal grouping term: g1/g2 \nGrouping terms may contain only variable names",
                     "combined by the interaction symbol ':'\n"))
})

test_that("Test that extract effects finds all response variables", {
 expect_equal(extract.effects(y1~x)$response, "y1")
 expect_equal(extract.effects(cbind(y1,y2)~x)$response, c("y1", "y2")) 
 expect_equal(extract.effects(cbind(y1,y2,y2)~x)$response, c("y1", "y2")) 
 expect_equal(extract.effects(y1+y2+y3~x)$response, c("y1", "y2", "y3")) 
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

test_that("Test that extract.time returns all desired variables", {
  expect_equal(extract.time(~1), list(time = "", group = "", all = ~1))
  expect_equal(extract.time(~tt), list(time = "tt", group = "", all = ~1 + tt)) 
  expect_equal(extract.time(~1|trait), list(time = "", group = "trait", all = ~1+trait)) 
  expect_equal(extract.time(~time|trait), 
               list(time = "time", group = "trait", all = ~1+time+trait)) 
  expect_equal(extract.time(~time|Site:trait),
               list(time = "time", group = "Site__trait", all = ~1+time+Site+trait))
  expect_error(extract.time(~t1+t2|g1), 
               "Autocorrelation structures may only contain 1 time variable")
  expect_error(extract.time(~1|g1/g2), 
    paste("Illegal grouping term: g1/g2 \nGrouping terms may contain only variable names",
          "combined by the interaction symbol ':'\n"))
})

test_that("Test that brm.update.formula returns correct formulas", {
  expect_equal(brm.update.formula(y~x, addition = list(se = ~I(sei+2))), y | se(I(sei+2)) ~ x)
  expect_equal(brm.update.formula(y~x, addition = list(se = ~sei, cens = ~censored)), 
               y | se(sei) | cens(censored) ~ x)
  expect_equal(brm.update.formula(y~x+z, partial = ~ a + I(a^2)), y ~ x+z+partial(a + I(a^2)))
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