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

test_that("Test that extract_effects finds all response variables", {
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
  expect_warning(update_formula(y~x, addition = list(se = ~I(sei+2))))
  expect_warning(update_formula(y~x, addition = list(se = ~sei, cens = ~censored)))
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

test_that("Test that check_family returns correct links", {
  expect_equal(check_family("gaussian")$link, "identity")
  expect_equal(check_family("weibull")$link, "log")
  expect_equal(check_family(binomial)$link, "logit")
  expect_equal(check_family(binomial("probit"))$link, "probit")
  expect_equal(check_family(c("acat", "cloglog"))$link, "cloglog")
})

test_that("Test that check_family return an error on wrong links", {
  expect_error(check_family(gaussian("logit")), 
               "logit is not a supported link for family gaussian")
  expect_error(check_family(poisson("inverse")), 
               "inverse is not a supported link for family poisson")
  expect_error(check_family(c("weibull", "sqrt")), 
               "sqrt is not a supported link for family weibull")
  expect_error(check_family(c("categorical","probit")), 
               "probit is not a supported link for family categorical")
})

test_that("Test that get_prior finds all classes for which priors can be specified", {
  expect_equal(get_prior(count ~ log_Base4_c * Trt_c + (1|patient) + (1+Trt_c|visit),
                       data = epilepsy, family = "poisson")$class,
               c(rep("b", 5), c("cor", "cor"), rep("sd", 6)))
  expect_equal(get_prior(rating ~ treat + period, partial = ~ carry, data = inhaler, 
                        family = "sratio", threshold = "equidistant")$class,
               c(rep("b", 5), "delta"))
})

test_that("Test that update_prior produces correct prior_frames", {
  prior <- list(b = "p1", sd = "p2", cor = "p3", b_Intercept = "p4",
                cor_visit = "p5", sd_visit_x = "p6", sd_visit = "p7", 
                sigma = "p8")
  result <- prior_frame(prior = paste0("p",1:8), 
                        class = c("b", "sd", "cor", "b", "cor", "sd", "sd", "sigma"),
                        coef = c(rep("", 3), "Intercept", "", "x", "", ""),
                        group = c(rep("", 4), rep("visit", 3), ""))
  expect_equal(update_prior(prior), result)
})

test_that("Test that gather_ranef works correctly", {
  data <- data.frame(g = 1:10, x = 11:20)
  target <- list(g = c("Intercept", "x"))
  attr(target$g, "levels") <- paste(1:10)
  expect_equal(gather_ranef(list(random = list(~1+x), group = "g"), data = data),
               target)
  expect_equal(gather_ranef(list()), list())
})
