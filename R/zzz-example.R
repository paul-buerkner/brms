# brmsfit_example <- brm(count ~ log_Age_c + (1|visit),
#                        data = data.frame(count = rpois(236, lambda = 20),
#                                          visit = rep(c(1:4), each = 59),
#                                          log_Age_c = rnorm(236)), 
#                        family = poisson("log"),
#                        prior = c(set_prior("normal(0,5)", class = "b"),
#                                  set_prior("cauchy(0,2)", class = "sd")),
#                        sample.prior = TRUE, n.warmup = 10, n.iter = 50)