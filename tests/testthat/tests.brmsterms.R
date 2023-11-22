context("Tests for formula parsing functions")

test_that("brmsterms finds all variables in very long formulas", {
  expect_equal(
    all.vars(brmsterms(t2_brand_recall ~ psi_expsi + psi_api_probsolv +
                                 psi_api_ident + psi_api_intere + psi_api_groupint)$all),
    all.vars(t2_brand_recall ~ t2_brand_recall + psi_expsi + psi_api_probsolv + psi_api_ident +
                 psi_api_intere + psi_api_groupint)
  )
})

test_that("brmsterms handles very long RE terms", {
  # tests issue #100
  covariate_vector <- paste0("xxxxx", 1:80, collapse = "+")
  formula <- paste(sprintf("y ~ 0 + trait + trait:(%s)", covariate_vector),
                   sprintf("(1+%s|id)", covariate_vector), sep = " + ")
  bterms <- brmsterms(as.formula(formula))
  expect_equal(bterms$dpars$mu$re$group, "id")
})

test_that("brmsterms correctly handles auxiliary parameter 'mu'", {
  bterms1 <- brmsterms(y ~ x + (x|g))
  bterms2 <- brmsterms(bf(y ~ 1, mu ~ x + (x|g)))
  expect_equal(bterms1$dpars$mu, bterms2$dpars$mu)

  # commented out for now as updating is not yet enabled
  # bterms1 <- brmsterms(bf(y ~ z + x + (x|g)))
  # bterms2 <- brmsterms(bf(y ~ z, lf(mu ~ x + (x|g))))
  # expect_equal(bterms1$dpars$mu, bterms2$dpars$mu)
  #
  # bterms1 <- brmsterms(bf(y ~ z, lf(mu ~ x + (x|g), cmc = FALSE)))
  # expect_true(!attr(bterms1$dpars$mu$fe, "cmc"))
  #
  # expect_error(brmsterms(bf(y ~ z, mu ~ x + (x|g), nl = TRUE)),
  #              "Cannot combine non-linear formulas")
})

test_that("brmsterms correctly check fixed auxiliary parameters", {
  bform <- bf(y~1, sigma = 4, family = gaussian)
  expect_true(is.brmsterms(brmsterms(bform)))
  bform <- bf(y~1, zi = 0.5, family = zero_inflated_beta())
  expect_true(is.brmsterms(brmsterms(bform)))
  bform <- bf(y~1, shape = -2, family = Gamma())
  expect_error(brmsterms(bform), "Parameter 'shape' must be positive")
  bform <- bf(y~1, quantile = 1.5, family = asym_laplace())
  expect_error(brmsterms(bform), "Parameter 'quantile' must be between 0 and 1")
})

test_that("check_re_formula returns correct REs", {
  old_form <- y ~ x + (1|patient) + (Trt_c|visit)
  form <- check_re_formula(~ (1 | visit), old_form)
  expect_equivalent(form, ~ (1 | gr(visit)))
  form <- check_re_formula(~ (1 + Trt_c|visit), old_form)
  expect_equivalent(form, ~ (1 + Trt_c | gr(visit)))
  form <- check_re_formula(~ (0 + Trt_c | visit) + (1|patient), old_form)
  expect_equivalent(form, ~ (1|gr(patient)) + (0 + Trt_c | gr(visit)))

  # checks for fix of issue #844
  old_form <- y ~ 0 + x1 + x2 + (0 + x1 + x2 | x3)
  expect_error(
    check_re_formula(~ (0 + x2 + x1 | x3), old_form),
    "Order of terms in 're_formula' should match the original order"
  )
})

test_that("update_re_terms works correctly", {
  expect_equivalent(update_re_terms(y ~ x, ~ (1|visit)), y ~ x)
  expect_equivalent(update_re_terms(y ~ x*z + (1+Trt_c|patient), ~ (1|patient)),
                    y ~ x*z + (1|gr(patient)))
  expect_equivalent(update_re_terms(y ~ x + (1|patient), ~ 1), y ~ x)
  expect_equivalent(update_re_terms(y ~ 1|patient, ~ 1), y ~ 1)
  expect_equivalent(update_re_terms(y ~ -1 + x + (1+visit|patient), NA),
                    y ~ -1 + x)
  expect_equivalent(update_re_terms(y ~ x + (1+visit|patient), NULL),
                    y ~ x + (1+visit|patient))
  expect_equivalent(update_re_terms(y ~ (1|patient), NA), y ~ 1)
  expect_equivalent(update_re_terms(y ~ x + (1+x|visit), ~ (1|visit)),
                    y ~ x + (1|gr(visit)))
  expect_equivalent(update_re_terms(y ~ x + (1|visit), ~ (1|visit) + (x|visit)),
                    y ~ x + (1|gr(visit)))
  expect_equal(update_re_terms(bf(y ~ x, sigma = ~ x + (x|g)), ~ (1|g)),
               bf(y ~ x, sigma = ~ x + (1|gr(g))))
  expect_equal(update_re_terms(bf(y ~ x, x ~ z + (1|g), nl = TRUE), ~ (1|g)),
               bf(y ~ x, x ~ z + (1|gr(g)), nl = TRUE))
})

test_that("unused variables are correctly incorporated", {
  bterms <- brmsterms(bf(y ~ 1, unused = ~ x))
  expect_true("x" %in% all.vars(bterms$allvars))
})
