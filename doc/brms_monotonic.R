params <-
  list(EVAL = TRUE)

## ---- SETTINGS-knitr, include=FALSE-----------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "jpeg",
  dpi = 100,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)
library(brms)
ggplot2::theme_set(theme_default())

## ---------------------------------------------------------------------------------------
income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(income_options, 100, TRUE),
  levels = income_options, ordered = TRUE
)
mean_ls <- c(30, 60, 70, 75)
ls <- mean_ls[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)

## ---- results='hide'--------------------------------------------------------------------
fit1 <- brm(ls ~ mo(income), data = dat)

## ---------------------------------------------------------------------------------------
summary(fit1)
plot(fit1, variable = "simo", regex = TRUE)
plot(conditional_effects(fit1))

## ---- results='hide'--------------------------------------------------------------------
dat$income_num <- as.numeric(dat$income)
fit2 <- brm(ls ~ income_num, data = dat)

## ---------------------------------------------------------------------------------------
summary(fit2)

## ---- results='hide'--------------------------------------------------------------------
contrasts(dat$income) <- contr.treatment(4)
fit3 <- brm(ls ~ income, data = dat)

## ---------------------------------------------------------------------------------------
summary(fit3)

## ---------------------------------------------------------------------------------------
loo(fit1, fit2, fit3)

## ---- results='hide'--------------------------------------------------------------------
prior4 <- prior(dirichlet(c(2, 1, 1)), class = "simo", coef = "moincome1")
fit4 <- brm(ls ~ mo(income),
  data = dat,
  prior = prior4, sample_prior = TRUE
)

## ---------------------------------------------------------------------------------------
summary(fit4)

## ---------------------------------------------------------------------------------------
plot(fit4, variable = "prior_simo", regex = TRUE, N = 3)

## ---------------------------------------------------------------------------------------
dat$age <- rnorm(100, mean = 40, sd = 10)

## ---- results='hide'--------------------------------------------------------------------
fit5 <- brm(ls ~ mo(income) * age, data = dat)

## ---------------------------------------------------------------------------------------
summary(fit5)
conditional_effects(fit5, "income:age")

## ---------------------------------------------------------------------------------------
dat$city <- rep(1:10, each = 10)
var_city <- rnorm(10, sd = 10)
dat$ls <- dat$ls + var_city[dat$city]

## ---- results='hide'--------------------------------------------------------------------
fit6 <- brm(ls ~ mo(income) * age + (mo(income) | city), data = dat)

## ---------------------------------------------------------------------------------------
summary(fit6)
