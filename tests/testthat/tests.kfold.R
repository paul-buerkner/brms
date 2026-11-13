context("Tests for kfold")

test_that("kfold warns about intercept priors with centered predictors", {
  user_intercept_prior <- prior(normal(0, 1), class = "Intercept")
  check_warning <- function(formula) {
    brms:::warn_kfold_intercept_prior(
      list(formula = formula, prior = user_intercept_prior)
    )
  }
  warning <- paste(
    "You are applying `kfold` to a model with a user-specified prior on",
    "the intercept and predictors or smooth terms that are centered or",
    "constrained using the supplied data.",
    "The meaning of this prior may vary across folds."
  )

  expect_warning(check_warning(bf(y ~ x)), warning, fixed = TRUE)
  expect_warning(
    check_warning(bf(mvbind(y1, y2) ~ x)), warning, fixed = TRUE
  )
  expect_warning(check_warning(bf(y ~ 0 + x)), NA)
  expect_warning(check_warning(bf(y ~ 0 + Intercept + x)), NA)
  expect_warning(check_warning(bf(y ~ x, center = FALSE)), NA)
  expect_warning(check_warning(bf(y ~ s(x))), warning, fixed = TRUE)
  expect_warning(check_warning(bf(y ~ 0 + Intercept + s(x))), NA)
  expect_warning(check_warning(bf(y ~ s(x), center = FALSE)), NA)

  expect_warning(check_warning(bf(y ~ 1, sigma ~ x)), warning, fixed = TRUE)
  expect_warning(
    check_warning(bf(y ~ 0 + x, sigma ~ z)), warning, fixed = TRUE
  )
  expect_warning(check_warning(bf(y ~ 1, sigma ~ 0 + x)), NA)
  expect_warning(
    check_warning(bf(y ~ 1, sigma ~ 0 + Intercept + x)), NA
  )
  expect_warning(
    check_warning(bf(y ~ 1, lf(sigma ~ x, center = FALSE))), NA
  )
  expect_warning(check_warning(bf(y ~ 1, sigma ~ s(x))), warning, fixed = TRUE)
  expect_warning(
    check_warning(bf(y ~ 0 + x, sigma ~ s(z))), warning, fixed = TRUE
  )
  expect_warning(check_warning(bf(y ~ 1, sigma ~ 0 + s(x))), NA)
  expect_warning(
    check_warning(bf(y ~ 1, lf(sigma ~ s(x), center = FALSE))), NA
  )
})
