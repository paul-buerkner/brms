set.seed(1234)
dat <- data.frame(
  count = rpois(236, lambda = 20),
  visit = rep(1:4, each = 59),
  patient = factor(rep(1:59, 4)),
  Age = rnorm(236), 
  Trt = factor(sample(0:1, 236, TRUE)),
  AgeSD = abs(rnorm(236, 1)),
  Exp = sample(1:5, 236, TRUE),
  volume = rnorm(236),
  gender = factor(c(rep("m", 30), rep("f", 29)))
)

dat2 <- data.frame(
  rating = sample(1:4, 50, TRUE), 
  subject = rep(1:10, 5),
  x1 = rnorm(50), 
  x2 = rnorm(50),
  x3 = rnorm(50)
)

library(brms)
brmsfit_example1 <- brm(
  bf(count ~ Trt*Age + mo(Exp) + s(Age) +
      offset(Age) + (1+Trt|visit), sigma ~ Trt),
  data = dat, family = student(), 
  autocor = cor_arma(~visit|patient, 1, 1),
  prior = c(set_prior("normal(0,5)", class = "b"),
            set_prior("cauchy(0,2)", class = "sd"),
            set_prior("normal(0,3)", dpar = "sigma")),
  sample_prior = TRUE, warmup = 150, iter = 200, chains = 2,
  save_dso = FALSE, testmode = TRUE
)

brmsfit_example2 <- brm(
  bf(count | weights(Exp) ~ inv_logit(a) * exp(b * Trt),
     a + b ~ Age + (1|ID1|patient), nl = TRUE),
  data = dat, family = Gamma(), 
  prior = c(set_prior("normal(2,2)", nlpar = "a"),
            set_prior("normal(0,3)", nlpar = "b")),
  sample_prior = TRUE, warmup = 150, iter = 200, chains = 2,
  save_dso = FALSE, testmode = TRUE
)

brmsfit_example3 <- brm(
  count ~ Trt*me(Age, AgeSD) + (1 + mmc(Age, volume) | mm(patient, visit)),
  data = dat[1:30, ], prior = prior(normal(0, 10)),
  warmup = 150, iter = 200, chains = 2, 
  save_mevars = TRUE, save_dso = FALSE, testmode = TRUE
)

brmsfit_example4 <- brm(
  bf(rating ~ x1 + cs(x2) + (cs(x2)||subject), disc ~ 1),
  data = dat2, family = sratio(),
  warmup = 150, iter = 200, chains = 2,
  save_dso = FALSE, testmode = TRUE
)

brmsfit_example5 <- brm(
  bf(count ~ Age + (1|gr(patient, by = gender)), mu2 ~ Age), 
  data = dat, family = mixture(gaussian, exponential),
  prior = c(prior(normal(0, 10), Intercept, dpar = mu1),
            prior(normal(0, 1), Intercept, dpar = mu2),
            prior(normal(0, 1), dpar = mu2)),
  warmup = 150, iter = 200, chains = 2,
  save_dso = FALSE, testmode = TRUE
)

brmsfit_example6 <- brm(
  bf(volume ~ Trt + gp(Age, by = Trt), family = gaussian()) +
    bf(count ~ Trt + Age, family = poisson()), 
  data = dat[1:40, ],
  prior = prior(normal(0, 0.25), lscale, resp = volume) +
    prior(normal(0, 10), sdgp, resp = volume),
  warmup = 150, iter = 200, chains = 2,
  save_dso = FALSE, testmode = TRUE
)

devtools::use_data(
  brmsfit_example1, brmsfit_example2, brmsfit_example3, 
  brmsfit_example4, brmsfit_example5, brmsfit_example6,
  internal = TRUE, overwrite = TRUE
)
