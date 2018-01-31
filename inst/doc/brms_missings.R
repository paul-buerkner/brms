params <-
structure(list(EVAL = TRUE), .Names = "EVAL")

## ---- SETTINGS-knitr, include=FALSE-----------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
  )

## ---------------------------------------------------------------------------------------
data("nhanes", package = "mice")
head(nhanes)

## ---------------------------------------------------------------------------------------
library(mice)
imp <- mice(nhanes, m = 5, print = FALSE)

## ---- results = 'hide', message = FALSE-------------------------------------------------
library(brms)
fit_imp1 <- brm_multiple(bmi ~ age*chl, data = imp, chains = 2)

## ---------------------------------------------------------------------------------------
summary(fit_imp1)

## ---------------------------------------------------------------------------------------
plot(fit_imp1, pars = "^b")

## ---------------------------------------------------------------------------------------
round(fit_imp1$rhats, 2)

## ---------------------------------------------------------------------------------------
marginal_effects(fit_imp1, "age:chl")

## ---- results = 'hide', message = FALSE-------------------------------------------------
bform <- bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE)
fit_imp2 <- brm(bform, data = nhanes)

## ---------------------------------------------------------------------------------------
summary(fit_imp2)
marginal_effects(fit_imp2, "age:chl", resp = "bmi")

