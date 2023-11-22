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
group <- rep(c("treat", "placebo"), each = 30)
symptom_post <- c(rnorm(30, mean = 1, sd = 2), rnorm(30, mean = 0, sd = 1))
dat1 <- data.frame(group, symptom_post)
head(dat1)

## ---- results='hide'--------------------------------------------------------------------
fit1 <- brm(bf(symptom_post ~ group, sigma ~ group),
            data = dat1, family = gaussian())

## ---- results='hide'--------------------------------------------------------------------
summary(fit1)
plot(fit1, N = 2, ask = FALSE)
plot(conditional_effects(fit1), points = TRUE)

## ---------------------------------------------------------------------------------------
hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_grouptreat) = 0")
hypothesis(fit1, hyp)

## ---------------------------------------------------------------------------------------
hyp <- "exp(sigma_Intercept + sigma_grouptreat) > exp(sigma_Intercept)"
(hyp <- hypothesis(fit1, hyp))
plot(hyp, chars = NULL)

## ---------------------------------------------------------------------------------------
zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")
head(zinb)

## ---- results='hide'--------------------------------------------------------------------
fit_zinb1 <- brm(count ~ persons + child + camper,
                 data = zinb, family = zero_inflated_poisson())

## ---------------------------------------------------------------------------------------
summary(fit_zinb1)
plot(conditional_effects(fit_zinb1), ask = FALSE)

## ---- results='hide'--------------------------------------------------------------------
fit_zinb2 <- brm(bf(count ~ persons + child + camper, zi ~ child),
                 data = zinb, family = zero_inflated_poisson())

## ---------------------------------------------------------------------------------------
summary(fit_zinb2)
plot(conditional_effects(fit_zinb2), ask = FALSE)

## ---------------------------------------------------------------------------------------
dat_smooth <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
head(dat_smooth[, 1:6])

## ---- results='hide'--------------------------------------------------------------------
fit_smooth1 <- brm(
  bf(y ~ s(x1) + s(x2) + (1|fac), sigma ~ s(x0) + (1|fac)),
  data = dat_smooth, family = gaussian(),
  chains = 2, control = list(adapt_delta = 0.95)
)

## ---------------------------------------------------------------------------------------
summary(fit_smooth1)
plot(conditional_effects(fit_smooth1), points = TRUE, ask = FALSE)

