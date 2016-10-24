test_that("check_prior performs correct renaming", {
  # ordering bug in devtools::check() for upper case letters
  prior <- check_prior(c(set_prior("lkj(0.5)", class = "cor", group = "subject"),
                         set_prior("normal(0,1)", class = "b")),
                       formula = bf(rating ~ treat + (0 + treat +  carry | subject)), 
                       data = inhaler)
  target <- prior_frame(prior = c("normal(0,1)", "lkj_corr_cholesky(0.5)"),
                        class = c("b", "L"), group = c("", "subject"))
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2L)
  
  prior <- check_prior(c(set_prior("lkj(5)", class = "rescor"),
                         set_prior("normal(0,1)", class = "b", 
                                   coef = "carry", nlpar = "rating")),
                       formula = bf(cbind(treat, rating) ~ 0 + carry), 
                       data = inhaler)
  target <- prior_frame(prior = c("lkj_corr_cholesky(5)", "normal(0,1)"),
                        class = c("Lrescor", "b"), coef = c("", "carry"),
                        nlpar = c("", "rating")) 
  expect_true(length(which(duplicated(rbind(prior, target)))) == 2L)
  
  expect_equivalent(check_prior(set_prior("normal(0,1)", class = "Intercept"),
                                formula = bf(rating ~ carry), data = inhaler, 
                                family = cumulative())[3, ],
                    prior_frame("normal(0,1)", class = "temp_Intercept"))
  
  expect_equivalent(check_prior(set_prior("normal(0,1)", class = "Intercept"),
                                formula = bf(rating ~ carry), data = inhaler, 
                                family = cumulative(),
                                threshold = "equidistant")[4, ],
                    prior_frame("normal(0,1)", class = "temp_Intercept"))
  
  expect_equivalent(check_prior(set_prior("normal(0,2)", "b", coef = "Intercept"),
                                formula = bf(rating ~ carry), data = inhaler, 
                                family = student())[5, ],
                    prior_frame("normal(0,2)", class = "temp_Intercept"))
})

test_that("check_prior accepts correct prior names", {
  cp <- check_prior(c(set_prior("normal(0,1)", coef = "carry"),
                      set_prior("cauchy(1,1)", coef = "treat")),
                    formula = bf(rating ~ -1 + treat + carry), data = inhaler)
  expect_equivalent(cp[2:3, ], prior_frame(c("normal(0,1)", "cauchy(1,1)"), 
                                 class = "b", coef = c("carry", "treat")))
  
  cp <- check_prior(c(set_prior("p1", class = "sd", coef = "sexfemale", 
                                group = "patient"),
                      set_prior("p2", class = "sd", coef = "age", 
                                group = "patient")),
                    formula = bf(time ~ age + (sex+age|patient)),  
                    family = exponential(), data = kidney)
  expect_equivalent(cp[9, ], prior_frame(prior = "p1", class = "sd", 
                                coef = "sexfemale", group = "patient")[1, ])
  
  expect_equivalent(check_prior(set_prior("cauchy(0,1)", class = "sigma"), 
                                formula = bf(rating ~ 1), family = gaussian(), 
                                data = inhaler)[2, ],
                    prior_frame("cauchy(0,1)", class = "sigma"))
  
  expect_equivalent(check_prior(c(set_prior("p1", class = "ar"),
                                  set_prior("p2", class = "ma")),
                                formula = bf(count ~ Trt_c), data = epilepsy, 
                                autocor = cor.arma(p = 1, q = 2))[c(1, 4), ],
                    prior_frame(c("p1", "p2"), class = c("ar", "ma"),
                                bound = "<lower=-1,upper=1>"))
  
  expect_equivalent(check_prior(set_prior("cauchy(0,1)", class = "sigmaLL"),
                                formula = bf(count ~ Trt_c), data = epilepsy,
                                autocor = cor_bsts())[4, ],
                    prior_frame("cauchy(0,1)", class = "sigmaLL"))
})

test_that("check_prior rejects incorrect prior names", {
  expect_message(check_prior(c(set_prior("p1", class = "Intercept"),
                               set_prior("p2", class = "b", coef = "age")),
                             family = acat(), data = inhaler,
                             formula = bf(rating ~ treat + (1+treat|subject))))
  expect_message(check_prior(c(set_prior("p1", class = "Intercept"),
                               set_prior("", class = "sd", group = "patient")),
                             formula = bf(rating ~ treat + (1+treat|subject)), 
                             family = student(), data = inhaler))
  expect_message(check_prior(set_prior("normal(0,1)", class = "ar"), 
                             formula = bf(count ~ log_Base4_c * Trt_c 
                             + (1+Trt_c|patient)), data = epilepsy))
})

test_that("check_prior returns increment_log_prob(.) whithout checking", {
  expect_equivalent(check_prior(c(set_prior("increment_log_prob(p1)"),
                                  set_prior("p2", class = "b")),
                                formula = bf(count ~ Trt_c), 
                                data = epilepsy)[c(1, 5), ],
                    prior_frame(c("p2", "increment_log_prob(p1)"), 
                                class = c("b", "")))
})

test_that("check_prior correctly validates priors for random effects", {
  expect_message(check_prior(set_prior("normal(0,1)", class = "sd", group = "g"),
                             formula = bf(count ~ (1|visit)), data = epilepsy),
                 "The following priors don't correspond")
  cp <- check_prior(set_prior("cauchy(0,1)", class = "sd", group = "visit"),
                    formula = bf(count ~ Trt_c + (1|visit)), 
                    data = epilepsy)
  expect_equal(cp$prior[4], "cauchy(0,1)")
  cp <- check_prior(set_prior("cauchy(0,1)", class = "sd", group = "visit"),
                    formula = bf(count ~ 1 + (1|visit) + (0+Trt_c|visit)), 
                    data = epilepsy)
  expect_equal(cp$prior[3], "cauchy(0,1)")
})

test_that("check_prior correctly validates priors for category specific effects", {
  prior <- c(set_prior("normal(0,1)", class = "b", coef = "carry"),
             set_prior("cauchy(1,1)", class = "b", coef = "treat"))
  cp <- check_prior(prior, formula = bf(rating ~ 1 + cse(treat + carry)), 
                    data = inhaler, family = cratio())
  target <- prior_frame(prior = c("normal(0,1)", "cauchy(1,1)"),
                        class = "b", coef = c("carry", "treat"))
  expect_equivalent(cp[2:3, ], target)
})

test_that("check_prior correctly validates priors for monotonic effects", {
  data <- data.frame(y = rpois(100, 10), x = rep(1:4, 25))
  prior <- c(set_prior("normal(0,1)", class = "b", coef = "x"),
             set_prior("dirichlet(c(1,0.5,2))", class = "simplex", coef = "x"))
  cp <- check_prior(prior, formula = bf(y ~ monotonic(x)), data = data,
                    family = poisson())
  target <- prior_frame(prior = c("normal(0,1)", "dirichlet(c(1,0.5,2))"),
                        class = c("b", "simplex"), coef = c("x", "x"))
  expect_equivalent(cp[2:3, ], target)
  expect_error(check_prior(set_prior("beta(1,1)", class = "simplex", coef = "x"), 
                           formula = bf(y ~ monotonic(x)), data = data),
               "'dirichlet' is the only valid prior for simplex parameters")
})

test_that("handle_special_priors handles horseshoe prior correctly", {
  prior <- set_prior("horseshoe(5)")
  prior <- handle_special_priors(c(prior))
  expect_equal(attr(prior, "hs_df"), 5)
  expect_equal(prior$prior[1], "normal(0, hs_local * hs_global)")
  expect_error(handle_special_priors(c(set_prior("horseshoe(-1)"))),
               "Degrees of freedom of the local priors")
  expect_error(handle_special_priors(c(set_prior("horseshoe(1, -1)"))),
               "Scale of the global prior")
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
  expect_equal(sort(unique(gp$nlpar)), c("", "a", "b"))
})

test_that("get_prior returns correct fixed effect names for GAMMs", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), 
                    z = rnorm(10), g = rep(1:2, 5))
  prior <- get_prior(y ~ z + s(x) + (1|g), data = dat)
  expect_equal(prior[prior$class == "b", ]$coef, 
               c("", "Intercept", "sx_1", "z"))
  prior <- get_prior(y ~ lp, nonlinear = lp ~ z + s(x) + (1|g), data = dat)
  expect_equal(prior[prior$class == "b", ]$coef, 
               c("", "Intercept", "sx_1", "z"))
})

test_that("get_prior returns correct prior names for auxiliary parameters", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), 
                    z = rnorm(10), g = rep(1:2, 5))
  prior <- get_prior(bf(y ~ 1, phi ~ z + (1|g)), data = dat, family = Beta())
  prior <- prior[prior$nlpar == "phi", ]
  pdata <- data.frame(class = rep(c("b", "sd"), each = 3), 
                      coef = c("", "Intercept", "z", "", "", "Intercept"),
                      group = c(rep("", 4), "g", "g"),
                      stringsAsFactors = FALSE)
  expect_equivalent(prior[, c("class", "coef", "group")], pdata)
})

test_that("get_prior returns global priors in multivariate models", {
  dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), 
                    x = rnorm(10), g = rep(1:2, 5))
  # MV normal
  prior <- get_prior(cbind(y1, y2) ~ x + (x|ID1|g), 
                     data = dat, family = gaussian())
  expect_equal(prior[prior$nlpar == "" & prior$class == "b", "coef"],
               c("", "Intercept", "x"))
  expect_equal(prior[prior$nlpar == "" & prior$class == "sd", "prior"],
               c("student_t(3, 0, 10)"))
  # categorical
  prior <- get_prior(y2 ~ x + (x|ID1|g), 
                     data = dat, family = categorical())
  expect_equal(prior[prior$nlpar == "" & prior$class == "b", "coef"],
               c("", "Intercept", "x"))
  expect_equal(prior[prior$nlpar == "" & prior$class == "sd", "prior"],
               c("student_t(3, 0, 10)"))
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

test_that("set_prior alias functions produce equivalent results", {
  expect_equal(set_prior("normal(0, 1)", class = "sd"),
               prior(normal(0, 1), class = sd))
  expect_equal(set_prior("normal(0, 1)", class = "sd", nlpar = "a"),
               prior(normal(0, 1), class = "sd", nlpar = a))
  expect_equal(set_prior("normal(0, 1)", class = "sd"),
               prior_string("normal(0, 1)", class = "sd"))
})
