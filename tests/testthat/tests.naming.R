test_that("Test that rename returns an error on duplicated names", {
  expect_error(rename(c(letters[1:4],"a()","a["), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: a, a(), a[")
  expect_error(rename(c("aDb","a/b","b"), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: aDb, a/b")
  expect_error(rename(c("log(a,b)","logab","bac","ba"), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: log(a,b), logab")
})

test_that("Test that check_prior performs correct renaming", {
  expect_equal(check_prior(list(cor = "normal(0,5)"), family = "student",
                           formula = rating ~ treat + (1+treat|subject), data = inhaler),
               list(L = "normal(0,5)"))
  expect_equal(check_prior(list(cor_subject = "lkj(0.5)", b = "normal(0,1)"), 
                           rating ~ treat + (1+treat|subject), data = inhaler),
               list(L_subject = "lkj_corr_cholesky(0.5)", b = "normal(0,1)"))
  expect_equal(check_prior(list(rescor = "lkj(5)", b_carry = "normal(0,1)"), 
                           cbind(rating, treat) ~ carry, data = inhaler),
               list(Lrescor = "lkj_corr_cholesky(5)", b_carry = "normal(0,1)"))
  expect_equal(check_prior(list(b_Intercept = "normal(0,1)"), family = "cumulative",
                           rating ~ carry, data = inhaler, threshold = "equidistant"),
               list(b_Intercept1 = "normal(0,1)"))
})

test_that("Test that check_prior is backwards compatible", {
  expect_equal(check_prior(list(L = "lkj_corr_cholesky(0.5)"), family = "student",
                           formula = rating ~ treat + (1+treat|subject), data = inhaler),
               list(L = "lkj_corr_cholesky(0.5)"))
  expect_equal(check_prior(list(Lrescor = "normal(0,5)", b_carry = "normal(0,1)"), 
                           cbind(rating, treat) ~ carry, data = inhaler),
               list(Lrescor = "normal(0,5)", b_carry = "normal(0,1)"))
})

test_that("Test that check_prior accepts correct prior names", {
  expect_equal(check_prior(list(b_Intercept = "normal(0,1)",  b_treat = "gamma(1,1)"), 
                 family = "acat", formula = rating ~ treat + (1+treat|subject), data = inhaler),
               list(b_Intercept = "normal(0,1)",  b_treat = "gamma(1,1)"))
  expect_equal(check_prior(list(b = "normal(0,1)",  sd = "gamma(1,1)"), 
                family = "exponential", formula = time ~ age + (1+age|patient), data = kidney),
               list(b = "normal(0,1)",  sd = "gamma(1,1)"))
  expect_equal(check_prior(list(sd_Intercept = "normal(0,1)",  sd_age = "gamma(1,1)"), 
                           family = "exponential", formula = time ~ age + (1+age|patient), data = kidney),
               list(sd_Intercept = "normal(0,1)",  sd_age = "gamma(1,1)"))
  expect_equal(check_prior(list(sigma = "cauchy(0,1)"), formula = rating ~ 1, 
                           family = "cauchy", data = inhaler),
               list(sigma = "cauchy(0,1)"))
  expect_equal(check_prior(list(nu = "cauchy(0,1)", b = "uniform(0,1)"), 
                 formula = rating ~ carry + (1|subject), family = "student", data = inhaler),
               list(nu = "cauchy(0,1)", b = "uniform(0,1)"))
  expect_equal(check_prior(list(shape = "cauchy(0,1)", b = "uniform(0,1)"), 
                           formula = time ~ sex + age, family = "gamma", data = kidney),
               list(shape = "cauchy(0,1)", b = "uniform(0,1)"))
  expect_equal(check_prior(list(delta = "gamma(1,1)"), formula = rating ~ 1, 
               family = "sratio", data = inhaler, threshold = "equidistant"),
               list(delta = "gamma(1,1)"))
  expect_equal(check_prior(list(ar = "normal(0,2)", ma = "student_t(1,2,3)"),
               formula = count ~ Trt_c, data = epilepsy, autocor = cor.arma(p = 1, q = 2)),
               list(ar = "normal(0,2)", ma = "student_t(1,2,3)"))
})

test_that("Test that check_priro rejects incorrect prior names", {
  expect_warning(check_prior(list(b_Intercept = "normal(0,1)",  b_age = "gamma(1,1)"), 
                 family = "acat", formula = rating ~ treat + (1+treat|subject), data = inhaler))
  expect_warning(check_prior(list(b_Intercept = "normal(0,1)",  sd_patient = "gamma(1,1)"), 
                 family = "cauchy", formula = rating ~ treat + (1+treat|subject), data = inhaler))
  expect_warning(check_prior(list(ar = "normal(0,1)",  b = "gamma(1,1)"), 
                 formula = count ~ log_Base4_c * Trt_c + (1+Trt_c|patient), data = epilepsy))
})
