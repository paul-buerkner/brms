test_that("check_prior performs correct renaming", {
  # ordering bug in devtools::check() for upper case letters
  prior <- check_prior(set_prior("normal(0,5)", class = "cor"),
                       formula = rating ~ 0 + treat + (0 + treat + carry | subject), 
                       data = inhaler, family = "student")
  target <- prior_frame(prior = "normal(0,5)", class = "L")
  expect_true(length(which(duplicated(rbind(prior, target)))) == 1)
  
  prior <- check_prior(c(set_prior("lkj(0.5)", class = "cor", group = "subject"),
                         set_prior("normal(0,1)", class = "b")),
                       formula = rating ~ treat + (0 + treat +  carry | subject), 
                       data = inhaler)
  target <- prior_frame(prior = c("normal(0,1)", "lkj_corr_cholesky(0.5)"),
                        class = c("b", "L"), group = c("", "1"))
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2)
  
  prior <- check_prior(c(set_prior("lkj(5)", class = "rescor"),
                         set_prior("normal(0,1)", class = "b", coef = "carry")),
                       formula = cbind(treat, rating) ~ 0 + carry, data = inhaler)
  target <- prior_frame(prior = c("lkj_corr_cholesky(5)", "normal(0,1)"),
                        class = c("Lrescor", "b"), coef = c("", "carry")) 
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2)
  
  expect_equivalent(check_prior(set_prior("normal(0,1)", class = "b", coef = "Intercept"),
                                formula = rating ~ carry, data = inhaler, 
                                family = "cumulative", link = "logit")[3, ],
                    prior_frame("normal(0,1)", class = "b_Intercept"))
  
  expect_equivalent(check_prior(set_prior("normal(0,1)", class = "b", coef = "Intercept"),
                                formula = rating ~ carry, data = inhaler, 
                                family = "cumulative", link = "logit",
                                threshold = "equidistant")[3, ],
                    prior_frame("normal(0,1)", class = "b_Intercept1"))
})


test_that("check_prior is backwards compatible", { 
  prior <- suppressWarnings(check_prior(
    list(b_carry = "normal(0,1)", nu = "gamma(1,1)"), 
    family = "student", formula = rating ~ carry + (1+treat|subject), 
    data = inhaler))
  target <- prior_frame(prior = c("normal(0,1)", "gamma(1,1)"),
                        class = c("b", "nu"), coef = c("carry", ""))
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2)
  
  prior <- suppressWarnings(check_prior(
    list(sd_subject_treat = "normal(0,1)", Lrescor = "lkj_corr_cholesky(1)"), 
    formula = cbind(rating, carry) ~ treat + (1+treat|subject), 
    data = inhaler, family = "gaussian"))
  target <- prior_frame(prior = c("normal(0,1)", "lkj_corr_cholesky(1)"),
                        class = c("sd", "Lrescor"), coef = c("treat", ""),
                        group = c("1", ""))
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2)
})

test_that("check_prior accepts correct prior names", {
  expect_equivalent(check_prior(c(set_prior("normal(0,1)", class = "b", coef = "carry"),
                                  set_prior("gamma(1,1)", class = "b", coef = "treat")),
                                formula = rating ~ -1 + treat + carry, data = inhaler)[c(2,3), ],
                    prior_frame(c("normal(0,1)", "gamma(1,1)"), class = "b",
                                coef = c("carry", "treat")))
  
  expect_equivalent(check_prior(c(set_prior("p1", class = "sd", coef = "sexfemale", group = "patient"),
                                  set_prior("p2", class = "sd", coef = "age", group = "patient")),
                                formula = time ~ age + (sex+age|patient),  
                                family = "exponential", data = kidney)[11, ],
                    prior_frame(prior = "p1", class = "sd", 
                                coef = "sexfemale", group = "1")[1, ])
  
  expect_equivalent(check_prior(set_prior("cauchy(0,1)", class = "sigma"), 
                                formula = rating ~ 1, family = "cauchy", data = inhaler)[4, ],
                    prior_frame("cauchy(0,1)", class = "sigma"))
  
  expect_equivalent(check_prior(c(set_prior("p1", class = "ar"),
                                  set_prior("p2", class = "ma")),
                                formula = count ~ Trt_c, data = epilepsy, 
                                autocor = cor.arma(p = 1, q = 2))[c(1, 6), ],
                    prior_frame(c("p1", "p2"), class = c("ar", "ma")))
})

test_that("check_prior rejects incorrect prior names", {
  expect_message(check_prior(c(set_prior("p1", class = "b", coef = "Intercept"),
                               set_prior("p2", class = "b", coef = "age")),
                             family = "acat", link = "logit", data = inhaler,
                             formula = rating ~ treat + (1+treat|subject)))
  expect_message(check_prior(c(set_prior("p1", class = "b", coef = "Intercept"),
                               set_prior("", class = "sd", group = "patient")),
                             formula = rating ~ treat + (1+treat|subject), 
                             family = "cauchy", data = inhaler))
  expect_message(check_prior(set_prior("normal(0,1)", class = "ar"), 
                             formula = count ~ log_Base4_c * Trt_c + (1+Trt_c|patient), 
                             data = epilepsy))
})

test_that("check_prior returns increment_log_prob(.) whithout checking", {
  expect_equivalent(check_prior(c(set_prior("increment_log_prob(p1)"),
                                  set_prior("p2", class = "b")),
                                formula = count ~ Trt_c, 
                                data = epilepsy)[c(1,7), ],
                    prior_frame(c("p2", "increment_log_prob(p1)"), 
                                class = c("b", "")))
}) 

test_that("handle_special_priors handles horseshoe prior correctly", {
  prior <- set_prior("horseshoe(5)")
  temp <- handle_special_priors(c(prior))
  expect_equal(temp$attrib$hs_df, 5)
  expect_equal(temp$prior$prior[1], "normal(0, hs_local * hs_global)")
  expect_error(handle_special_priors(c(prior, set_prior("dist()", coef = "a"))))
  expect_error(handle_special_priors(c(set_prior("horseshoe(b5)"))),
               "degrees of freedom of horseshoe prior must be a positive number")
})

test_that("get_prior finds all classes for which priors can be specified", {
  expect_equal(get_prior(count ~ log_Base4_c * Trt_c + (1|patient) + (1+Trt_c|visit),
                         data = epilepsy, family = "poisson")$class,
               c(rep("b", 5), c("cor", "cor"), rep("sd", 6)))
  expect_equal(get_prior(rating ~ treat + period, partial = ~ carry, data = inhaler, 
                         family = "sratio", threshold = "equidistant")$class,
               c(rep("b", 5), "delta"))
})

test_that("update_prior produces correct prior_frames", {
  prior <- list(b = "p1", sd = "p2", cor = "p3", b_Intercept = "p4",
                cor_visit = "p5", sd_visit_x = "p6", sd_visit = "p7", 
                sigma = "p8")
  result <- prior_frame(prior = paste0("p",1:8), 
                        class = c("b", "sd", "cor", "b", "cor", "sd", "sd", "sigma"),
                        coef = c(rep("", 3), "Intercept", "", "x", "", ""),
                        group = c(rep("", 4), rep("visit", 3), ""))
  expect_equal(update_prior(prior), result)
})
