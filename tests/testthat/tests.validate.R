test_that("Test that extract_effects finds all variables in very long formulas", {
  expect_equal(extract_effects(t2_brand_recall ~ psi_expsi + psi_api_probsolv + 
                                 psi_api_ident + psi_api_intere + psi_api_groupint)$all, 
               t2_brand_recall ~ psi_expsi + psi_api_probsolv + psi_api_ident + psi_api_intere + psi_api_groupint)
})

test_that("Test that extract_effects finds all random effects and grouping factors", {
  expect_equal(extract_effects(y ~ a + (1+x|g1) + x + (1|g2) + z)$random, list(~1 + x, ~1))
  expect_equal(extract_effects(y ~ (1+x|g1) + x + (1|g2))$group, c("g1", "g2"))
  expect_equal(extract_effects(y ~ (1+x|g1:g2) + x + (1|g1))$group, c("g1", "g1:g2"))
  expect_error(extract_effects(y ~ (1+x|g1/g2) + x + (1|g1)))
})

test_that("Test that extract_effects accepts || syntax", {
  expect_equal(extract_effects(y ~ a + (1+x||g1) + (1+z|g2))$cor, c(FALSE,TRUE))
  expect_equal(extract_effects(y ~ a + (1+x||g2))$random, list(~1 + x))
  expect_equal(extract_effects(y ~ (1+x||g1) + x + (1||g2))$group, c("g1", "g2"))
  expect_equal(extract_effects(y ~ (1+x||g1:g2))$group, c("g1:g2"))
  expect_error(extract_effects(y ~ (1+x||g1/g2) + x + (1|g1)))
})

test_that("Test that extract effects finds all response variables", {
  expect_equal(extract_effects(y1~x)$response, "y1")
  expect_equal(extract_effects(cbind(y1,y2)~x)$response, c("y1", "y2")) 
  expect_equal(extract_effects(cbind(y1,y2,y2)~x)$response, c("y1", "y2")) 
  expect_equal(extract_effects(y1+y2+y3~x)$response, c("y1", "y2", "y3")) 
})

test_that("Test that extract_effects handles addition arguments correctly", {
  expect_equal(extract_effects(y | se(I(a+2)) ~ x, family = "gaussian")$se, ~ .se(I(a+2)))
  expect_equal(extract_effects(y | se(I(a+2)) ~ x, family = "gaussian")$all, y ~ x + a)
  expect_equal(extract_effects(y | weights(1/n) ~ x, family = "gaussian")$weights, ~ .weights(1/n))
  expect_equal(extract_effects(y | se(a+2) | cens(log(b)) ~ x, family = "gaussian")$cens, ~ .cens(log(b)))
  expect_equal(extract_effects(y | trials(10) ~ x, family = "binomial")$trials, 10)
  expect_equal(extract_effects(y | cat(cate) ~ x, family = "cumulative")$cat, ~ .cat(cate))
  expect_equal(extract_effects(y | cens(cens^2) ~ z, family = "weibull")$cens, ~ .cens(cens^2))
  expect_equal(extract_effects(y | cens(cens^2) ~ z + (x|patient), family = "weibull")$all, 
               y ~ z + x + patient + cens)
})

test_that("Test that extract_time returns all desired variables", {
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
               paste("Illegal grouping term: g1/g2 \n",
                     "may contain only variable names combined by the symbol ':'\n"))
})

test_that("Test that update_formula returns correct formulas", {
  expect_equal(update_formula(y~x, addition = list(se = ~I(sei+2))), y | se(I(sei+2)) ~ x)
  expect_equal(update_formula(y~x, addition = list(se = ~sei, cens = ~censored)), 
               y | se(sei) | cens(censored) ~ x)
  expect_equal(update_formula(y~x+z, partial = ~ a + I(a^2)), y ~ x+z+partial(a + I(a^2)))
})

test_that("Test that get_group_formula rejects incorrect grouping terms", {
  expect_error(get_group_formula("|g1/g2"), 
               paste("Illegal grouping term: g1/g2 \n",
                     "may contain only variable names combined by the symbol ':'"))
  expect_error(get_group_formula("||g1/g2"), 
               paste("Illegal grouping term: g1/g2 \n",
                     "may contain only variable names combined by the symbol ':'"))
})

test_that("Test that check_prior performs correct renaming", {
  expect_equal(check_prior(set_prior("normal(0,5)", class = "cor"),
                           formula = rating ~ treat + (1+treat|subject), 
                           data = inhaler, family = "student"),
               prior_frame(prior = "normal(0,5)", class = "L"))
  
  expect_equal(check_prior(list(set_prior("lkj(0.5)", class = "cor", group = "subject"),
                                set_prior("normal(0,1)", class = "b")),
                           rating ~ treat + (1+treat|subject), data = inhaler),
               prior_frame(prior = c("lkj_corr_cholesky(0.5)", "normal(0,1)"),
                          class = c("L", "b"), group = c("1", "")))
  
  expect_equal(check_prior(list(set_prior("lkj(5)", class = "rescor"),
                                set_prior("normal(0,1)", class = "b", coef = "carry")),
                           cbind(rating, treat) ~ carry, data = inhaler),
               prior_frame(prior = c("lkj_corr_cholesky(5)", "normal(0,1)"),
                           class = c("Lrescor", "b"), coef = c("", "carry")))
  
  expect_equal(check_prior(set_prior("normal(0,1)", class = "b", coef = "Intercept"),
                           formula = rating ~ carry, data = inhaler, family = "cumulative"),
               prior_frame("normal(0,1)", class = "b_Intercept"))
  
  expect_equal(check_prior(set_prior("normal(0,1)", class = "b", coef = "Intercept"),
                           formula = rating ~ carry, data = inhaler, 
                           family = "cumulative", threshold = "equidistant"),
               prior_frame("normal(0,1)", class = "b_Intercept1"))
})

# TODO 
#test_that("Test that check_prior is backwards compatible", { 
#  expect_equal(check_prior(list(L = "lkj_corr_cholesky(0.5)"), family = "student",
#                           formula = rating ~ treat + (1+treat|subject), data = inhaler),
#               list(L = "lkj_corr_cholesky(0.5)"))
#  expect_equal(check_prior(list(Lrescor = "normal(0,5)", b_carry = "normal(0,1)"), 
#                           cbind(rating, treat) ~ carry, data = inhaler),
#               list(Lrescor = "normal(0,5)", b_carry = "normal(0,1)"))
#})

test_that("Test that check_prior accepts correct prior names", {
  expect_equal(check_prior(list(set_prior("normal(0,1)", class = "b", coef = "Intercept"),
                                set_prior("gamma(1,1)", class = "b", coef = "treat")),
                           data = inhaler, formula = rating ~ treat + (1+treat|subject)),
               prior_frame(c("normal(0,1)", "gamma(1,1)"), class = "b",
                           coef = c("Intercept", "treat")))

  expect_equal(check_prior(list(set_prior("p1", class = "sd", coef = "Intercept",
                                          group = "patient"),
                                set_prior("p2", class = "sd", coef = "age",
                                          group = "patient")),
                           family = "exponential", data = kidney,
                           formula = time ~ age + (1+age|patient)),
               prior_frame(prior = c("p1", "p2"), class = "sd", 
                           coef = c("Intercept", "age"), group = "1"))
  
  expect_equal(check_prior(set_prior("cauchy(0,1)", class = "sigma"), 
                           formula = rating ~ 1, family = "cauchy", data = inhaler),
               prior_frame("cauchy(0,1)", class = "sigma"))
  
  expect_equal(check_prior(list(set_prior("p1", class = "ar"),
                                set_prior("p2", class = "ma")),
                           formula = count ~ Trt_c, data = epilepsy, 
                           autocor = cor.arma(p = 1, q = 2)),
               prior_frame(c("p1", "p2"), class = c("ar", "ma")))
})

test_that("Test that check_prior rejects incorrect prior names", {
  expect_message(check_prior(list(set_prior("p1", class = "b", coef = "Intercept"),
                                  set_prior("p2", class = "b", coef = "age")),
                             family = "acat", data = inhaler,
                             formula = rating ~ treat + (1+treat|subject)))
  expect_error(check_prior(list(set_prior("p1", class = "b", coef = "Intercept"),
                                  set_prior("", class = "sd", group = "patient")),
                             formula = rating ~ treat + (1+treat|subject), 
                             family = "cauchy", data = inhaler))
  expect_message(check_prior(set_prior("normal(0,1)", class = "ar"), 
                             formula = count ~ log_Base4_c * Trt_c + (1+Trt_c|patient), 
                             data = epilepsy))
})

test_that("Test that check_family rejects invalid families", {
  expect_error(check_family("multigaussian"),
               "family 'multigaussian' is deprecated. Use family 'gaussian' instead")
  expect_error(check_family("ordinal"),
               "ordinal is not a valid family")
})

test_that("Test that link4familys returns correct links", {
  expect_equal(link4family("gaussian"), "identity")
  expect_equal(link4family("weibull"), "log")
  expect_equal(link4family("binomial"), "logit")
  expect_equal(link4family(c("binomial", "probit")), "probit")
  expect_equal(link4family(c("acat", "cloglog")), "cloglog")
  expect_warning(link4family(c("poisson", "sqrt")), "poisson model with sqrt link may not be uniquely identified")
})

test_that("Test that link4familys return an error on wrong links", {
  expect_error(link4family(c("gaussian","logit")), "logit is not a valid link for family gaussian")
  expect_error(link4family(c("poisson","inverse")), "inverse is not a valid link for family poisson")
  expect_error(link4family(c("weibull","sqrt")), "sqrt is not a valid link for family weibull")
  expect_error(link4family(c("categorical","probit")), "probit is not a valid link for family categorical")
})

test_that("Test that parnames.formula finds all classes for which priors can be specified", {
  expect_equal(parnames(count ~ log_Base4_c * Trt_c + (1|patient) + (1+Trt_c|visit),
                       data = epilepsy, family = "poisson")$class,
               c(rep("b", 5), rep("sd", 6), c("cor", "cor")))
  expect_equal(parnames(rating ~ treat + period, partial = ~ carry, data = inhaler, 
                        family = "sratio", threshold = "equidistant")$class,
               c(rep("b", 5), "delta"))
})

test_that("Test that update_deprecated_prior produces correct prior_frames", {
  prior <- list(b = "p1", sd = "p2", cor = "p3", b_Intercept = "p4",
                cor_visit = "p5", sd_visit_z_x = "p6", sigma = "p7")
  result <- prior_frame(prior = paste0("p",1:7), 
                        class = c("b", "sd", "cor", "b", "cor", "sd", "sigma"),
                        coef = c(rep("", 3), "Intercept", "", "z_x", ""),
                        group = c(rep("", 4), rep("visit", 2), ""))
  expect_equal(update_deprecated_prior(prior),
               result)
})
