test_that("parse_bf finds all variables in very long formulas", {
  expect_equal(parse_bf(t2_brand_recall ~ psi_expsi + psi_api_probsolv + 
                                 psi_api_ident + psi_api_intere + psi_api_groupint)$all, 
               t2_brand_recall ~ t2_brand_recall + psi_expsi + psi_api_probsolv + psi_api_ident + 
                 psi_api_intere + psi_api_groupint)
})

test_that("parse_bf handles very long RE terms", {
  # tests issue #100
  covariate_vector <- paste0("xxxxx", 1:80, collapse = "+")
  formula <- paste(sprintf("y ~ 0 + trait + trait:(%s)", covariate_vector),
                   sprintf("(1+%s|id)", covariate_vector), sep = " + ")
  bterms <- parse_bf(as.formula(formula))
  expect_equal(bterms$dpars$mu$re$group, "id")
})

test_that("parse_bf correctly handles auxiliary parameter 'mu'", {
  bterms1 <- parse_bf(y ~ x + (x|g))
  bterms2 <- parse_bf(bf(y~1, mu ~ x + (x|g)))
  expect_equal(bterms1$dpars$mu, bterms2$dpars$mu)
  expect_error(parse_bf(bf(y ~ z, mu ~ x + (x|g))),
               "All 'mu' parameters are specified")
})

test_that("parse_bf correctly check fixed auxiliary parameters", {
  bform <- bf(y~1, sigma = 4, family = gaussian)
  expect_true(is.brmsterms(parse_bf(bform)))
  bform <- bf(y~1, zi = 0.5, family = zero_inflated_beta())
  expect_true(is.brmsterms(parse_bf(bform)))
  bform <- bf(y~1, shape = -2, family = Gamma())
  expect_error(parse_bf(bform), "Parameter 'shape' must be positive")
  bform <- bf(y~1, quantile = 1.5, family = asym_laplace())
  expect_error(parse_bf(bform), "Parameter 'quantile' must be between 0 and 1")
})

test_that("check_re_formula returns correct REs", {
  old_form <- y ~ x + (1|patient) + (Trt_c|visit)
  form <- check_re_formula(~ (1 | visit), old_form)
  expect_equivalent(form, ~ (1 | gr(visit)))
  form <- check_re_formula(~ (1 + Trt_c|visit), old_form)
  expect_equivalent(form, ~ (1 + Trt_c | gr(visit)))
  form <- check_re_formula(~ (0 + Trt_c | visit) + (1|patient), old_form)
  expect_equivalent(form, ~ (1|gr(patient)) + (0 + Trt_c | gr(visit)))
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

test_that("exclude_pars returns expected parameter names", {
  ranef <- data.frame(id = c(1, 1, 2), group = c("g1", "g1", "g2"),
                       gn = c(1, 1, 2), coef = c("x", "z", "x"), 
                       cn = c(1, 2, 1), nlpar = "", 
                       cor = c(TRUE, TRUE, FALSE))
  empty_effects <- structure(list(), class = "brmsterms")
  ep <- exclude_pars(empty_effects, ranef = ranef)
  expect_true(all(c("r_1", "r_2") %in% ep))
  ep <- exclude_pars(empty_effects, ranef = ranef, save_ranef = FALSE)
  expect_true("r_1_1" %in% ep)
  
  ranef$nlpar <- c("a", "a", "")
  ep <- exclude_pars(empty_effects, ranef = ranef, save_ranef = FALSE)
  expect_true(all(c("r_1_a_1", "r_1_a_2") %in% ep))
  bterms <- parse_bf(y ~ x + s(z))
  data <- data.frame(y = rnorm(20), x = rnorm(20), z = rnorm(20))
  expect_true("zs_1_1" %in% exclude_pars(bterms, data))
  bterms <- parse_bf(bf(y ~ eta, eta ~ x + s(z), family = gaussian(), nl = TRUE))
  expect_true("zs_eta_1_1" %in% exclude_pars(bterms, data))
})
