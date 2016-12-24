test_that("extract_effects finds all variables in very long formulas", {
  expect_equal(extract_effects(t2_brand_recall ~ psi_expsi + psi_api_probsolv + 
                                 psi_api_ident + psi_api_intere + psi_api_groupint)$all, 
               t2_brand_recall ~ t2_brand_recall + psi_expsi + psi_api_probsolv + psi_api_ident + 
                 psi_api_intere + psi_api_groupint)
})

test_that("extract_effects finds all random effects terms", {
  random <- extract_effects(y ~ a + (1+x|g1) + x + (1|g2) + z)$random
  expect_equal(random$form, list(~1 + x, ~1))
  expect_equal(random$group, c("g1", "g2"))
  expect_equal(extract_effects(y ~ (1+x|g1:g2) + (1|g1))$random$group, 
               c("g1", "g1:g2"))
  expect_equal(extract_effects(y ~ (1|g1))$random$group, c("g1"))
  expect_equal(extract_effects(y ~ log(x) + (1|g1))$random$group, c("g1"))
  expect_equal(extract_effects(y ~ (x+z):v + (1+x|g1))$random$form, list(~1+x))
  expect_equal(extract_effects(y ~ v + (1+x|g1/g2/g3))$random$group, 
               c("g1", "g1:g2", "g1:g2:g3"))
  expect_equal(extract_effects(y ~ v + (1+x|g1/g2) + (1|g3))$random$form, 
               list(~1+x, ~1+x, ~1))
  expect_equal(extract_effects(y ~ (1+x||g1/g2) + (1|g3))$random$cor, 
               c(FALSE, FALSE, TRUE))
  expect_equal(extract_effects(y ~ 1|g1)$random$group, "g1")
  expect_equal(extract_effects(y ~ x + (1+x|g1+g2))$random$group,
               c("g1", "g2"))
})

test_that("extract_effects finds all response variables", {
  expect_equal(extract_effects(y1~x)$response, "y1")
  expect_equal(extract_effects(cbind(y1,y2)~x)$response, 
               c("y1", "y2")) 
  expect_equal(extract_effects(cbind(y1,y2,y2)~x)$response, 
               c("y1", "y2", "y21")) 
  expect_equal(extract_effects(y1+y2+y3~x)$response, "y1") 
  expect_equal(extract_effects(y1/y2 ~ (1|g))$response, "y1")
  expect_equal(extract_effects(cbind(y1/y2,y2,y3*3) ~ (1|g))$response,
               c("response1", "y2", "response3"))
})

test_that("extract_effects also saves untransformed variables", {
  ee <- extract_effects(y ~ as.numeric(x) + (as.factor(z) | g))
  expect_equivalent(ee$allvars, 
                    y ~ y + as.numeric(x) + x + as.factor(z) + g + z)
})

test_that("extract_effects finds all spline terms", {
  ee <- extract_effects(y ~ s(x) + t2(z) + v)
  expect_equal(all.vars(ee$fixed), c("y", "v"))
  expect_equivalent(ee$gam, ~ s(x) + t2(z))
  ee <- extract_effects(bf(y ~ lp, lp ~ s(x) + t2(z) + v, nl = TRUE))
  expect_equal(all.vars(ee$nlpars[[1]]$fixed), "v")
  expect_equivalent(ee$nlpars[[1]]$gam, ~ s(x) + t2(z))
  expect_error(extract_effects(y ~ s(x) + te(z) + v), 
               "splines 'te' and 'ti' are not yet implemented")
})

test_that("extract_effects handles very long RE terms", {
  # tests issue #100
  covariate_vector <- paste0("xxxxx", 1:80, collapse = "+")
  formula <- paste(sprintf("y ~ 0 + trait + trait:(%s)", covariate_vector),
                   sprintf("(1+%s|id)", covariate_vector), sep = " + ")
  ee <- extract_effects(formula = as.formula(formula))
  expect_equal(ee$random$group, "id")
})

test_that("amend_formula returns correct formulas", {
  expect_warning(uf <- amend_formula(y ~ x + z, partial = ~ a + I(a^2)))
  expect_equal(uf, y ~ x + z + cs(a + I(a^2)))
})

test_that("check_re_formula returns correct REs", {
  old_form <- y ~ x + (1|patient) + (Trt_c|visit)
  form <- check_re_formula(~(1|visit), old_form)
  expect_equivalent(form, ~(1|visit))
  form <- check_re_formula(~(1+Trt_c|visit), old_form)
  expect_equivalent(form, ~(1+Trt_c|visit))
  form <- check_re_formula(~(0+Trt_c|visit) + (1|patient), old_form)
  expect_equivalent(form, ~ (1|patient) + (0+Trt_c | visit))
})

test_that("update_re_terms works correctly", {
  expect_equivalent(update_re_terms(y ~ x, ~ (1|visit)), y ~ x)
  expect_equivalent(update_re_terms(y ~ x + (1+Trt_c|patient), ~ (1|patient)), 
                    y ~ x + (1|patient))
  expect_equivalent(update_re_terms(y ~ x + (1|patient), ~ 1), 
                    y ~ x)
  expect_equivalent(update_re_terms(y ~ x + (1+visit|patient), NA), 
                    y ~ x)
  expect_equivalent(update_re_terms(y ~ x + (1+visit|patient), NULL), 
                    y ~ x + (1+visit|patient))
  expect_equivalent(update_re_terms(y ~ (1|patient), NA), y ~ 1)
  expect_equivalent(update_re_terms(y ~ x + (1+x|visit), ~ (1|visit)), 
                    y ~ x + (1|visit))
  expect_equivalent(update_re_terms(y ~ x + (1|visit), ~ (1|visit) + (x|visit)),
                    y ~ x + (1|visit))
  expect_equal(update_re_terms(bf(y ~ x, sigma = ~ x + (x|g)), ~ (1|g)),
               bf(y ~ x, sigma = ~ x + (1|g)))
  expect_equal(update_re_terms(bf(y ~ x, x ~ z + (1|g), nl = TRUE), ~ (1|g)),
               bf(y ~ x, x ~ z + (1|g), nl = TRUE))
})

test_that("amend_terms performs expected changes to terms objects", {
  expect_equal(amend_terms("a"), NULL)
  expect_equal(amend_terms(y~x), terms(y~x))
  form <- structure(y ~ x, rsv_intercept = TRUE)
  expect_true(attr(amend_terms(form), "rm_intercept"))
  form <- structure(y ~ 0 + main, forked = TRUE, old_mv = TRUE)
  expect_true(attr(amend_terms(form), "rm_intercept"))
  form <- structure(y ~ main, forked = TRUE, old_mv = TRUE)
  expect_error(amend_terms(form), "formula may not contain an intercept")
  form <- structure(y ~ main + trait, forked = TRUE, old_mv = TRUE)
  expect_error(amend_terms(form), "formula may not contain variable 'trait'")
})

test_that("exclude_pars returns expected parameter names", {
  ranef <- data.frame(id = c(1, 1, 2), group = c("g1", "g1", "g2"),
                       gn = c(1, 1, 2), coef = c("x", "z", "x"), 
                       cn = c(1, 2, 1), nlpar = "", 
                       cor = c(TRUE, TRUE, FALSE))
  ep <- exclude_pars(list(), ranef = ranef)
  expect_true(all(c("r_1", "r_2") %in% ep))
  ep <- exclude_pars(list(), ranef = ranef, save_ranef = FALSE)
  expect_true("r_1_1" %in% ep)
  
  ranef$nlpar <- c("a", "a", "")
  ep <- exclude_pars(list(), ranef = ranef, save_ranef = FALSE)
  expect_true(all(c("r_1_a_1", "r_1_a_2") %in% ep))
  effects <- extract_effects(y ~ x + s(z))
  data <- data.frame(y = rnorm(20), x = rnorm(20), z = rnorm(20))
  expect_true("zs_1_1" %in% exclude_pars(effects, data))
  effects <- extract_effects(bf(y ~ eta, eta ~ x + s(z), nl = TRUE))
  expect_true("zs_eta_1_1" %in% exclude_pars(effects, data))
})
