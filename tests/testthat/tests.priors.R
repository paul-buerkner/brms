test_that("check_prior performs correct renaming", {
  # ordering bug in devtools::check() for upper case letters
  prior <- check_prior(c(set_prior("lkj(0.5)", class = "cor", group = "subject"),
                         set_prior("normal(0,1)", class = "b")),
                       formula = rating ~ treat + (0 + treat +  carry | subject), 
                       data = inhaler)
  target <- prior_frame(prior = c("normal(0,1)", "lkj_corr_cholesky(0.5)"),
                        class = c("b", "L"), group = c("", "subject"))
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2)
  
  prior <- check_prior(c(set_prior("lkj(5)", class = "rescor"),
                         set_prior("normal(0,1)", class = "b", coef = "carry")),
                       formula = cbind(treat, rating) ~ 0 + carry, data = inhaler)
  target <- prior_frame(prior = c("lkj_corr_cholesky(5)", "normal(0,1)"),
                        class = c("Lrescor", "b"), coef = c("", "carry")) 
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2)
  
  expect_equivalent(check_prior(set_prior("normal(0,1)", class = "Intercept"),
                                formula = rating ~ carry, data = inhaler, 
                                family = cumulative())[3, ],
                    prior_frame("normal(0,1)", class = "temp_Intercept"))
  
  expect_equivalent(check_prior(set_prior("normal(0,1)", class = "Intercept"),
                                formula = rating ~ carry, data = inhaler, 
                                family = cumulative(),
                                threshold = "equidistant")[4, ],
                    prior_frame("normal(0,1)", class = "temp_Intercept1"))
})

test_that("check_prior accepts correct prior names", {
  cp <- check_prior(c(set_prior("normal(0,1)", coef = "carry"),
                      set_prior("cauchy(1,1)", coef = "treat")),
                    formula = rating ~ -1 + treat + carry, data = inhaler)
  expect_equivalent(cp[2:3, ], prior_frame(c("normal(0,1)", "cauchy(1,1)"), 
                                 class = "b", coef = c("carry", "treat")))
  
  cp <- check_prior(c(set_prior("p1", class = "sd", coef = "sexfemale", 
                                group = "patient"),
                      set_prior("p2", class = "sd", coef = "age", 
                                group = "patient")),
                    formula = time ~ age + (sex+age|patient),  
                    family = exponential(), data = kidney)
  expect_equivalent(cp[9, ], prior_frame(prior = "p1", class = "sd", 
                                coef = "sexfemale", group = "patient")[1, ])
  
  expect_equivalent(check_prior(set_prior("cauchy(0,1)", class = "sigma"), 
                                formula = rating ~ 1, family = cauchy(), 
                                data = inhaler)[2, ],
                    prior_frame("cauchy(0,1)", class = "sigma"))
  
  expect_equivalent(check_prior(c(set_prior("p1", class = "ar"),
                                  set_prior("p2", class = "ma")),
                                formula = count ~ Trt_c, data = epilepsy, 
                                autocor = cor.arma(p = 1, q = 2))[c(1, 4), ],
                    prior_frame(c("p1", "p2"), class = c("ar", "ma"),
                                bound = "<lower=-1,upper=1>"))
})

test_that("check_prior rejects incorrect prior names", {
  expect_message(check_prior(c(set_prior("p1", class = "Intercept"),
                               set_prior("p2", class = "b", coef = "age")),
                             family = acat(), data = inhaler,
                             formula = rating ~ treat + (1+treat|subject)))
  expect_message(check_prior(c(set_prior("p1", class = "Intercept"),
                               set_prior("", class = "sd", group = "patient")),
                             formula = rating ~ treat + (1+treat|subject), 
                             family = cauchy(), data = inhaler))
  expect_message(check_prior(set_prior("normal(0,1)", class = "ar"), 
                             formula = count ~ log_Base4_c * Trt_c 
                             + (1+Trt_c|patient), data = epilepsy))
})

test_that("check_prior returns increment_log_prob(.) whithout checking", {
  expect_equivalent(check_prior(c(set_prior("increment_log_prob(p1)"),
                                  set_prior("p2", class = "b")),
                                formula = count ~ Trt_c, 
                                data = epilepsy)[c(1, 6), ],
                    prior_frame(c("p2", "increment_log_prob(p1)"), 
                                class = c("b", "")))
})

test_that("check_prior correctly validates priors for random effects", {
  expect_message(check_prior(set_prior("normal(0,1)", class = "sd", group = "g"),
                             formula = count ~ (1|visit), data = epilepsy),
                 "The following priors don't correspond")
  cp <- check_prior(set_prior("cauchy(0,1)", class = "sd", group = "visit"),
                    formula = count ~ Trt_c + (1|visit), 
                    data = epilepsy)
  expect_equal(cp$prior[4], "cauchy(0,1)")
  cp <- check_prior(set_prior("cauchy(0,1)", class = "sd", group = "visit"),
                    formula = count ~ (1|visit) + (0+Trt_c|visit), 
                    data = epilepsy)
  expect_equal(cp$prior[3], "cauchy(0,1)")
})

test_that("check_prior correctly validates prior for category specific effects", {
  prior <- c(set_prior("normal(0,1)", class = "b", coef = "carry"),
             set_prior("cauchy(1,1)", class = "b", coef = "treat"))
  cp <- check_prior(prior, formula = rating ~ 1 + cse(treat + carry), 
                    data = inhaler, family = cratio())
  target <- prior_frame(prior = c("normal(0,1)", "cauchy(1,1)"),
                        class = "bp", coef = c("carry", "treat"))
  expect_equivalent(cp[3:4, ], target)
})

test_that("handle_special_priors handles horseshoe prior correctly", {
  prior <- set_prior("horseshoe(5)")
  temp <- handle_special_priors(c(prior))
  expect_equal(temp$attrib$hs_df, 5)
  expect_equal(temp$prior$prior[1], "normal(0, hs_local * hs_global)")
  expect_error(handle_special_priors(c(prior, set_prior("dist()", coef = "a"))))
  expect_error(handle_special_priors(c(set_prior("horseshoe(b5)"))),
               "degrees of freedom of horseshoe prior must be a positive")
})

test_that("get_prior finds all classes for which priors can be specified", {
  expect_equal(sort(get_prior(count ~ log_Base4_c * Trt_c 
                              + (1|patient) + (1+Trt_c|visit),
                              data = epilepsy, family = "poisson")$class),
               sort(c(rep("b", 5), c("cor", "cor"), "Intercept", rep("sd", 6))))
  expect_equal(sort(get_prior(rating ~ treat + period + cse(carry),
                         data = inhaler, family = "sratio", 
                         threshold = "equidistant")$class),
               sort(c(rep("b", 5), "delta", "Intercept")))
})

test_that("update_prior_frame works correctly", {
  prior <- prior_frame(group = c("", "1", "1", "2"))
  ranef <- list(g1 = "Intercept", g2 = c("Intercept", "x"))
  target <- prior_frame(group = c("", "g1", "g1", "g2"))
  expect_equal(update_prior_frame(prior, ranef), target)
  prior <- prior_frame(group = c("", "1", "1", "2", "3"))
  expect_warning(update_prior_frame(prior, ranef),
                 "Returning to default priors")
})

test_that("print for class brmsprior works correctly", {
  expect_output(print(set_prior("normal(0,1)")), fixed = TRUE,
                "b ~ normal(0,1)")
  expect_output(print(set_prior("normal(0,1)", coef = "x")), 
                "b_x ~ normal(0,1)", fixed = TRUE)
  expect_output(print(set_prior("cauchy(0,1)", class = "sd", group = "x")), 
                "sd_x ~ cauchy(0,1)", fixed = TRUE)
  expect_output(print(set_prior("increment_log_prob(normal_log(0,1))")), 
                "increment_log_prob(normal_log(0,1))", fixed = TRUE)
})

test_that("get_prior returns correct nlpar names for random effects pars", {
  # reported in issue #47
  data <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:2, 5))
  gp <- get_prior(y ~ a - b^x, data = data, nonlinear = a + b ~ (1+x|g))
  expect_equal(unique(gp$nlpar), c("", "a", "b"))
})

test_that("check_prior_content returns expected errors and warnings", {
  prior <- c(set_prior("", lb = 0), set_prior("gamma(0,1)", coef = "x"))
  expect_silent(check_prior_content(prior))
  prior <- c(set_prior("gamma(1,1)", class = "delta"))
  expect_silent(check_prior_content(prior, family = cumulative()))
  expect_warning(check_prior_content(prior, family = acat()),
                 "no natural lower bound")
  prior <- c(set_prior("uniform(0,5)", class = "sd"))
  expect_warning(check_prior_content(prior), "no natural upper bound")
  prior <- c(set_prior("normal(0,2)", class = "ar", lb = "0"))
  expect_warning(check_prior_content(prior), "autocorrelation parameters")
})
