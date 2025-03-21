set.seed(1234)
N <- 40
dat <- data.frame(
  count = rpois(N, lambda = 20),
  visit = factor(rep(1:4, each = N/4)),
  patient = factor(rep(1:(N/4), 4)),
  Age = rnorm(N),
  Trt = factor(sample(0:1, N, TRUE)),
  AgeSD = abs(rnorm(N, 1)),
  Exp = factor(sample(1:5, N, TRUE), ordered = TRUE),
  volume = rnorm(N),
  gender = factor(c(rep("m", N/8), rep("f", N/8)))
)

dat2 <- data.frame(
  rating = sample(1:4, N, TRUE),
  subject = rep(1:(N/5), 5),
  x1 = rnorm(N),
  x2 = rnorm(N),
  x3 = rnorm(N)
)

warmup <- 50
iter <- 75
chains <- 1
stan_model_args <- list(save_dso = FALSE)

library(brms)
brmsfit_example1 <- brm(
  bf(count ~ Trt*Age + mo(Exp) + s(Age) + volume +
      offset(Age) + (1+Trt|visit) + arma(visit, patient),
     sigma ~ Trt),
  data = dat, family = student(),
  prior = set_prior("normal(0,2)", class = "b") +
    set_prior("cauchy(0,2)", class = "sd") +
    set_prior("normal(0,3)", dpar = "sigma"),
  sample_prior = TRUE,
  warmup = warmup, iter = iter, chains = chains,
  stan_model_args = stan_model_args, rename = FALSE
)

brmsfit_example2 <- brm(
  bf(count | weights(AgeSD) ~ 1 / (1 + exp(-a)) * exp(b * Trt),
     a + b ~ Age + (1|ID1|patient), nl = TRUE),
  data = dat, family = Gamma("identity"),
  prior = set_prior("normal(2,2)", nlpar = "a") +
    set_prior("normal(0,3)", nlpar = "b"),
  sample_prior = TRUE,
  warmup = warmup, iter = iter, chains = chains,
  stan_model_args = stan_model_args, rename = FALSE
)

brmsfit_example3 <- brm(
  count ~ Trt*me(Age, AgeSD) + (1 + mmc(Age, volume) | mm(patient, visit)),
  data = dat[1:30, ], prior = prior(normal(0, 10)),
  save_mevars = TRUE,
  warmup = warmup, iter = iter, chains = chains,
  stan_model_args = stan_model_args, rename = FALSE
)

brmsfit_example4 <- brm(
  bf(rating ~ x1 + cs(x2) + (cs(x2)||subject), disc ~ 1),
  data = dat2, family = sratio(),
  warmup = warmup, iter = iter, chains = chains,
  stan_model_args = stan_model_args, rename = FALSE
)

brmsfit_example5 <- brm(
  bf(count ~ Age + (1|gr(patient, by = gender)), mu2 ~ Age),
  data = dat, family = mixture(gaussian, exponential),
  prior = prior(normal(0, 10), Intercept, dpar = mu1) +
    prior(normal(0, 1), Intercept, dpar = mu2) +
    prior(normal(0, 1), dpar = mu2),
  warmup = warmup, iter = iter, chains = chains,
  stan_model_args = stan_model_args, rename = FALSE
)

brmsfit_example6 <- brm(
  bf(volume ~ Trt + gp(Age, by = Trt, gr = TRUE), family = gaussian()) +
    bf(count ~ Trt + Age, family = poisson()),
  data = dat[1:40, ],
  prior = prior(normal(0, 0.25), lscale, resp = volume) +
    prior(normal(0, 10), sdgp, resp = volume),
  warmup = warmup, iter = iter, chains = chains,
  stan_model_args = stan_model_args, rename = FALSE
)

brmsfit_example7 <- SW(brm(
  formula = count ~ Trt + (1 | patient) + (1 + Trt | visit),
  data = dat[1:40, ],
  prior = c(prior(normal(0, 1), class = sd, tag = "prior_tag1"),
            prior(normal(0, 5), class = b, tag = "prior_tag2"),
            prior(normal(0, 0.5), coef = "Trt1", tag = "prior_tag3"),
            prior(normal(0, 10), class = "Intercept", tag = "prior_tag4"),
            prior(lkj_corr_cholesky(3), class = "L", group = "visit", tag = "prior_tag5")
            ),
  family = poisson()))


# easy loading of unchanged models to avoid refitting all of them
# brmsfit_example1 <- brms:::brmsfit_example1
# brmsfit_example2 <- brms:::brmsfit_example2
# brmsfit_example3 <- brms:::brmsfit_example3
# brmsfit_example4 <- brms:::brmsfit_example4
# brmsfit_example5 <- brms:::brmsfit_example5
# brmsfit_example6 <- brms:::brmsfit_example6
# brmsfit_example7 <- brms:::brmsfit_example7

usethis::use_data(
  brmsfit_example1, brmsfit_example2, brmsfit_example3,
  brmsfit_example4, brmsfit_example5, brmsfit_example6,
  brmsfit_example7,
  internal = TRUE, overwrite = TRUE
)
