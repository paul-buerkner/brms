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
  expect_error(extract_effects(y ~ (1+x|g1/g2) + x + (1|g1)))
  expect_error(extract_effects(y ~ 1|g1),
               "Random effects terms should be enclosed in brackets")
})

test_that("extract_effects accepts || syntax", {
  random <- extract_effects(y ~ a + (1+x||g1) + (1+z|g2))$random
  target <- data.frame(group = c("g1", "g2"), cor = c(FALSE, TRUE),
                       stringsAsFactors = FALSE)
  target$form <- list(~1+x, ~1+z)
  expect_equal(random, target)
  expect_equal(extract_effects(y ~ (1+x||g1:g2))$random$group, c("g1:g2"))
  expect_error(extract_effects(y ~ (1+x||g1/g2) + x + (1|g1)))
})

test_that("extract_effects finds all response variables", {
  expect_equal(extract_effects(y1~x)$response, "y1")
  expect_equal(extract_effects(cbind(y1,y2)~x)$response, 
               c("y1", "y2")) 
  expect_equal(extract_effects(cbind(y1,y2,y2)~x)$response, 
               c("y1", "y2", "y2")) 
  expect_equal(extract_effects(y1+y2+y3~x)$response, "y1") 
  expect_equal(extract_effects(y1/y2 ~ (1|g))$response, "y1")
  expect_equal(extract_effects(cbind(y1/y2,y2,y3*3) ~ (1|g))$response,
               c("response1", "y2", "response3"))
})

test_that("extract_effects handles addition arguments correctly", {
  expect_equal(extract_effects(y | se(I(a+2)) ~ x, family = gaussian())$se, 
               ~ .se(I(a+2)))
  expect_equal(extract_effects(y | se(I(a+2)) ~ x, family = gaussian())$all, 
               y ~ y + a + x)
  expect_equal(extract_effects(y | weights(1/n) ~ x, 
                               family = gaussian())$weights, 
               ~ .weights(1/n))
  expect_equal(extract_effects(y | se(a+2) | cens(log(b)) ~ x, 
                               family = gaussian())$cens, 
               ~ .cens(log(b)))
  expect_equal(extract_effects(y | se(a+2) + cens(log(b)) ~ x, 
                               family = gaussian())$se, 
               ~ .se(a+2))
  expect_equal(extract_effects(y | trials(10) ~ x, family = binomial())$trials, 
               10)
  expect_equal(extract_effects(y | cat(cate) ~ x, family = cumulative())$cat, 
               ~ .cat(cate))
  expect_equal(extract_effects(y | cens(cens^2) ~ z, family = weibull())$cens, 
               ~ .cens(cens^2))
  expect_equal(extract_effects(y | cens(cens^2) ~ z + (x|patient), 
                               family = weibull())$all, 
               y ~ y + cens + z + x + patient)
  expect_equal(extract_effects(resp | disp(a + b) ~ x, 
                               family = gaussian())$disp,
               ~.disp(a + b))
})

test_that("extract_effects accepts complicated random terms", {
  expect_equal(extract_effects(y ~ x + (I(as.numeric(x)-1) | z))$random$form,
               list(~I(as.numeric(x) - 1)))
  expect_equal(extract_effects(y ~ x + (I(exp(x)-1) + I(x/y) | z))$random$form,
               list(~I(exp(x)-1) + I(x/y)))
})

test_that("extract_effects accepts calls to the poly function", {
  expect_equal(extract_effects(y ~ z + poly(x, 3))$all,
               y ~ y + z + x + poly(x, 3))
})

test_that("extract_effects finds all variables in non-linear models", {
  nonlinear <- list(a ~ z1 + (1|g1), b ~ z2 + (z3|g2))
  ee <- extract_effects(y ~ a - b^x, nonlinear = nonlinear)
  expect_equal(ee$all, y ~ y + x + z1 + g1 + z2 + z3 + g2)
})

test_that("extract_effects rejects REs in non-linear formulas", {
  expect_error(extract_effects(y ~ exp(-x/a) + (1|g), nonlinear = a ~ 1),
               "Random effects in non-linear models", fixed = TRUE)
})

test_that("nonlinear_effects rejects invalid non-linear models", {
  expect_error(nonlinear_effects(list(a ~ 1, b ~ 1), model = y ~ a^x),
               "missing in formula: b")
  expect_error(nonlinear_effects(a~1, model = y ~ a^x),
               "Argument 'nonlinear' must be a list of formulas")
  expect_error(nonlinear_effects(list( ~ 1, a ~ 1), model = y ~ a),
               "Non-linear formulas must be two-sided")
  expect_error(nonlinear_effects(list(a + b ~ 1), model = y ~ exp(-x)),
               "LHS of non-linear formula must contain exactly one variable")
  expect_error(nonlinear_effects(list(a.b ~ 1), model = y ~ a^x),
               "not contain dots or underscores")
  expect_error(nonlinear_effects(list(a_b ~ 1), model = y ~ a^(x+b)),
               "not contain dots or underscores")
})

test_that("nonlinear_effect accepts valid non-linear models", {
  nle <- nonlinear_effects(list(a ~ 1 + (1+x|origin), b ~ 1 + z), y ~ b - a^x)
  expect_equal(names(nle), c("a", "b"))
  expect_equal(nle[["a"]]$all, ~x + origin)
  expect_equal(nle[["b"]]$all, ~z)
  expect_equal(nle[["a"]]$random$form[[1]], ~1+x)
})

test_that("nonlinear2list works correctly", {
  expect_equal(nonlinear2list(a ~ 1), list(a ~ 1))
  expect_equal(nonlinear2list(a + alpha ~ x + (x|g)), 
               list(a ~ x + (x|g), alpha ~ x + (x|g)))
  expect_equal(nonlinear2list(list(a ~ 1, b ~ 1 + z)),
               list(a ~ 1, b ~ 1 + z))
  expect_equal(nonlinear2list(NULL), NULL)
  expect_error(nonlinear2list(1), "invalid 'nonlinear' argument")
})

test_that("extract_time returns all desired variables", {
  expect_equal(extract_time(~1), list(time = "", group = "", all = ~1))
  expect_equal(extract_time(~tt), list(time = "tt", group = "", all = ~1 + tt)) 
  expect_equal(extract_time(~1|trait), list(time = "", group = "trait", all = ~1+trait)) 
  expect_equal(extract_time(~time|trait), 
               list(time = "time", group = "trait", all = ~1+time+trait)) 
  expect_equal(extract_time(~time|Site:trait),
               list(time = "time", group = "Site:trait", all = ~1+time+Site+trait))
  expect_error(extract_time(~t1+t2|g1), 
               "Autocorrelation structures may only contain 1 time variable")
  expect_error(extract_time(~1|g1/g2), 
               paste("Illegal grouping term: g1/g2"))
})

test_that("update_formula returns correct formulas", {
  expect_warning(uf <- update_formula(y ~ x + z, partial = ~ a + I(a^2)))
  expect_equal(uf, y ~ x + z + cse(a + I(a^2)))
})

test_that("get_effect works correctly", {
  effects <- extract_effects(y ~ a - b^x, 
               nonlinear = list(a ~ z, b ~ v + monotonous(z)))
  expect_equivalent(get_effect(effects), list(y ~ a - b^x, ~ z, ~ v))
  expect_equivalent(get_effect(effects, "mono"), list(NULL, NULL, ~ z))
  effects <- extract_effects(y ~ x + z + (1|g))
  expect_equivalent(get_effect(effects), list(y ~ x + z))
})

test_that("get_group_formula rejects incorrect grouping terms", {
  expect_error(get_group_formula("|g1/g2"), 
               paste("Illegal grouping term: g1/g2 \n",
                     "may contain only variable names combined by the symbol ':'"))
  expect_error(get_group_formula("||g1/g2"), 
               paste("Illegal grouping term: g1/g2 \n",
                     "may contain only variable names combined by the symbol ':'"))
})

test_that("check_re_formula returns correct REs", {
  old_ranef = list(patient = c("Intercept"), visit = c("Trt_c", "Intercept"))
  expect_equivalent(check_re_formula(~(1|visit), old_ranef = old_ranef, 
                                data = epilepsy),
                    list(visit = "Intercept"))
  expect_equivalent(check_re_formula(~(1+Trt_c|visit), old_ranef = old_ranef, 
                                data = epilepsy),
                    list(visit = c("Intercept", "Trt_c")))
  expect_equivalent(check_re_formula(~(0+Trt_c|visit) + (1|patient), 
                                     old_ranef = old_ranef, data = epilepsy),
                    list(patient = "Intercept", visit = "Trt_c"))
})

test_that("check_re_formula rejects invalid re_formulae", {
  old_ranef = list(patient = c("Intercept"), visit = c("Trt_c", "Intercept"))
  expect_error(check_re_formula(~ visit + (1|visit), old_ranef = old_ranef, 
                                data = epilepsy),
               "fixed effects are not allowed in re_formula")
  expect_error(check_re_formula(count ~ (1+Trt_c|visit), old_ranef = old_ranef, 
                                data = epilepsy),
               "re_formula must be one-sided")
  expect_error(check_re_formula(~(1|Trt_c), old_ranef = old_ranef, 
                                data = epilepsy),
               "Invalid grouping factors detected: Trt_c")
  expect_error(check_re_formula(~(1+Trt_c|patient), old_ranef = old_ranef, 
                                data = epilepsy),
               "Invalid random effects detected for grouping factor patient: Trt_c")
})

test_that("update_re_terms works correctly", {
  expect_equivalent(update_re_terms(y ~ x, ~ (1|visit)), y ~ x + (1|visit))
  expect_equivalent(update_re_terms(y ~ x + (1|patient), ~ (1|visit)), 
                    y ~ x + (1|visit))
  expect_equivalent(update_re_terms(y ~ x + (1|patient), ~ 1), 
                    y ~ x)
  expect_equivalent(update_re_terms(y ~ x + (1+visit|patient), NA), 
                    y ~ x)
  expect_equivalent(update_re_terms(y ~ x + (1+visit|patient), NULL), 
                    y ~ x + (1+visit|patient))
  expect_equivalent(update_re_terms(y ~ (1|patient), NA), y ~ 1)
  expect_equivalent(update_re_terms(y ~ x + (1+x|visit), ~ (1|visit)), 
                    y ~ x + (1|visit))
})

test_that("amend_terms performs expected changes to terms objects", {
  expect_equal(amend_terms("a"), NULL)
  expect_equal(amend_terms(y~x), terms(y~x))
  form <- structure(y~x, rsv_intercept = TRUE)
  expect_equal(attr(amend_terms(form), "rm_intercept"), TRUE)
  t <- amend_terms(y ~ 0 + main + main:x + spec + spec:z, forked = TRUE)
  expect_equal(attr(t, "intercept"), 1)
  expect_equal(attr(t, "rm_intercept"), TRUE)
  expect_error(amend_terms(y ~ main, forked = TRUE), "intercept")
  expect_error(amend_terms(y ~ 0 + main + trait, forked = TRUE), 
               "trait")
})

test_that("gather_ranef works correctly", {
  data <- data.frame(g = 1:10, x = 11:20, y = 1:10)
  target <- list(g = c("Intercept", "x"))
  attr(target$g, "levels") <- paste(1:10)
  attr(target$g, "group") <- "g"
  attr(target$g, "cor") <- FALSE
  expect_equal(gather_ranef(extract_effects(y~(1+x||g)), data = data),
               target)
  expect_equal(gather_ranef(list()), list())
})

test_that("check_brm_input returns correct warnings and errors", {
  expect_error(check_brm_input(list(chains = 3, cluster = 2)), 
               "chains must be a multiple of cluster", fixed = TRUE)
  x <- list(family = weibull(), inits = "random", chains = 1, cluster = 1,
            algorithm = "sampling")
  expect_warning(check_brm_input(x))
  x$family <- inverse.gaussian()
  expect_warning(check_brm_input(x))
  x$family <- poisson("sqrt")
  expect_warning(check_brm_input(x))
})

test_that("exclude_pars returns expected parameter names", {
  expect_true("eta_a" %in% exclude_pars(list(), nlpars = "a"))
  ranef <- list(g = structure(c("x", "z"), cor = TRUE))
  expect_true(all(c("r_1_1", "r_1_2") %in% exclude_pars(ranef)))
  expect_true("r_a_1_1" %in% exclude_pars(ranef, nlpars = "a"))
  expect_true("r_1" %in% exclude_pars(ranef, save_ranef = FALSE))
})
