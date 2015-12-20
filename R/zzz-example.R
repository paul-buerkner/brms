# Uncomment the code below to enable automated unit tests of S3 methods

brmsfit_example <- brm(count ~ Trt_c + offset(log_Age_c) + (1+Trt_c|visit),
                       data = data.frame(count = rpois(236, lambda = 20),
                                         visit = rep(c(1:4), each = 59),
                                         log_Age_c = rnorm(236),
                                         Trt_c = rnorm(236)), 
                       family = poisson("log"),
                       prior = c(set_prior("normal(0,5)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       sample.prior = TRUE, n.warmup = 10, n.iter = 50)
brmsfit_example$fit@stanmodel <- new("stanmodel")