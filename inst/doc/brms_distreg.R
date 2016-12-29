params <-
structure(list(EVAL = FALSE), .Names = "EVAL")

## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = params$EVAL,
  dev = "png",
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
  )

## ------------------------------------------------------------------------
#  group <- rep(c("treat", "placebo"), each = 30)
#  symptom_post <- c(rnorm(30, mean = 1, sd = 2), rnorm(30, mean = 0, sd = 1))
#  dat1 <- data.frame(group, symptom_post)
#  head(dat1)

## ---- results='hide'-----------------------------------------------------
#  fit1 <- brm(bf(symptom_post ~ group, sigma ~ group),
#              data = dat1, family = gaussian())

## ---- results='hide'-----------------------------------------------------
#  summary(fit1)
#  plot(fit1)
#  plot(marginal_effects(fit1), points = TRUE)

## ------------------------------------------------------------------------
#  hyp <- c("exp(sigma_Intercept) = 0",
#           "exp(sigma_Intercept + sigma_grouptreat) = 0")
#  hypothesis(fit1, hyp)

## ------------------------------------------------------------------------
#  hyp <- "exp(sigma_Intercept + sigma_grouptreat) > exp(sigma_Intercept)"
#  (hyp <- hypothesis(fit1, hyp))
#  plot(hyp, chars = NULL)

## ------------------------------------------------------------------------
#  zinb <- read.csv("http://www.ats.ucla.edu/stat/data/fish.csv")
#  head(zinb)

## ---- results='hide'-----------------------------------------------------
#  fit_zinb1 <- brm(count ~ persons + child + camper,
#                   data = zinb, family = zero_inflated_poisson())

## ------------------------------------------------------------------------
#  summary(fit_zinb1)
#  plot(fit_zinb1)
#  plot(marginal_effects(fit_zinb1), ask = FALSE)

## ---- results='hide'-----------------------------------------------------
#  fit_zinb2 <- brm(bf(count ~ persons + child + camper, zi ~ child),
#                   data = zinb, family = zero_inflated_poisson())

## ------------------------------------------------------------------------
#  summary(fit_zinb2)
#  plot(fit_zinb2, N = 6)
#  plot(marginal_effects(fit_zinb2), ask = FALSE)

## ------------------------------------------------------------------------
#  library(mgcv)
#  dat_smooth <- gamSim(eg = 6, n = 200, scale = 2)
#  head(dat_smooth)

## ---- results='hide'-----------------------------------------------------
#  fit_smooth1 <- brm(bf(y ~ s(x1) + s(x2) + (1|fac), sigma ~ s(x0) + (1|fac)),
#                     data = dat_smooth, family = gaussian(),
#                     chains = 2, control = list(adapt_delta = 0.95))

## ------------------------------------------------------------------------
#  summary(fit_smooth1)
#  plot(fit_smooth1, ask = FALSE)
#  plot(marginal_effects(fit_smooth1), points = TRUE, ask = FALSE)

