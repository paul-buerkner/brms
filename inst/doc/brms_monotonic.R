params <-
structure(list(EVAL = FALSE), .Names = "EVAL")

## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = FALSE, # params$EVAL,
  dev = "png",
  dpi = 150,
  fig.asp = 0.618,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
  )

## ------------------------------------------------------------------------
#  income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
#  income <- factor(sample(income_options, 100, TRUE),
#                   levels = income_options, ordered = TRUE)
#  mean_ls <- c(30, 60, 70, 75)
#  ls <- mean_ls[income] + rnorm(100, sd = 7)
#  dat <- data.frame(income, ls)

## ---- results='hide', message=FALSE, warning = FALSE---------------------
#  library(brms)
#  fit1 <- brm(ls ~ monotonic(income), data = dat)

## ------------------------------------------------------------------------
#  summary(fit1)
#  plot(fit1, pars = "simplex")
#  plot(marginal_effects(fit1))

## ---- results='hide', message=FALSE, warning = FALSE---------------------
#  dat$income_num <- as.numeric(dat$income)
#  fit2 <- brm(ls ~ income_num, data = dat)

## ------------------------------------------------------------------------
#  summary(fit2)

## ---- results='hide', message=FALSE, warning = FALSE---------------------
#  contrasts(dat$income) <- contr.treatment(4)
#  fit3 <- brm(ls ~ income, data = dat)

## ------------------------------------------------------------------------
#  summary(fit3)

## ------------------------------------------------------------------------
#  LOO(fit1, fit2, fit3)

## ---- results='hide', message=FALSE, warning = FALSE---------------------
#  prior4 <- prior(dirichlet(c(2, 1, 1)), class = "simplex", coef = "income")
#  fit4 <- brm(ls ~ monotonic(income), data = dat,
#             prior = prior4, sample_prior = TRUE)

## ------------------------------------------------------------------------
#  summary(fit4)

## ------------------------------------------------------------------------
#  plot(fit4, pars = "prior_simplex", N = 3)

## ------------------------------------------------------------------------
#  dat$city <- rep(1:10, each = 10)
#  var_city <- rnorm(10, sd = 10)
#  dat$ls <- dat$ls + var_city[dat$city]

## ---- results='hide', message=FALSE, warning = FALSE---------------------
#  fit5 <- brm(ls ~ mo(income) + (1 | city) + (mo(income) | city), data = dat)

## ------------------------------------------------------------------------
#  summary(fit5)

