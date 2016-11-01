## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
#  income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
#  income <- factor(sample(income_options, 100, TRUE),
#                   levels = income_options, ordered = TRUE)
#  mean_ls <- c(30, 60, 70, 75)
#  ls <- mean_ls[income] + rnorm(100, sd = 7)
#  dat <- data.frame(income, ls)

## ------------------------------------------------------------------------
#  library(brms)
#  fit1 <- brm(ls ~ monotonic(income), data = dat)

## ------------------------------------------------------------------------
#  summary(fit1)
#  plot(fit1, N = 6)
#  plot(marginal_effects(fit1), points = TRUE)

## ------------------------------------------------------------------------
#  dat$income_num <- as.numeric(dat$income)
#  fit2 <- brm(ls ~ income_num, data = dat)
#  summary(fit2)

## ------------------------------------------------------------------------
#  contrasts(dat$income) <- contr.treatment(4)
#  fit3 <- brm(ls ~ income, data = dat)
#  summary(fit3)

## ------------------------------------------------------------------------
#  LOO(fit1, fit2, fit3)

## ------------------------------------------------------------------------
#  prior4 <- prior(dirichlet(c(2, 1, 1)), class = "simplex", coef = "income")
#  fit4 <- brm(ls ~ monotonic(income), data = dat,
#              prior = prior4, sample_prior = TRUE)
#  summary(fit4)

## ------------------------------------------------------------------------
#  plot(fit4, pars = "simplex", N = 6)

## ------------------------------------------------------------------------
#  dat$city <- rep(1:10, each = 10)
#  var_city <- rnorm(10, sd = 10)
#  dat$ls <- dat$ls + var_city[dat$city]

## ------------------------------------------------------------------------
#  fit5 <- brm(ls ~ mono(income) + (1 | city) + (mono(income) | city), data = dat)

## ------------------------------------------------------------------------
#  summary(fit5)
#  plot(fit5, ask = FALSE)

