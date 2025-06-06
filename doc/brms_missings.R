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
data("nhanes", package = "mice")
head(nhanes)

## ---------------------------------------------------------------------------------------
library(mice)
imp <- mice(nhanes, m = 5, print = FALSE)

## ---- results = 'hide', message = FALSE-------------------------------------------------
fit_imp1 <- brm_multiple(bmi ~ age * chl, data = imp, chains = 2)

## ---------------------------------------------------------------------------------------
summary(fit_imp1)

## ---------------------------------------------------------------------------------------
plot(fit_imp1, variable = "^b", regex = TRUE)

## ---------------------------------------------------------------------------------------
round(fit_imp1$rhats, 2)

## ---------------------------------------------------------------------------------------
conditional_effects(fit_imp1, "age:chl")

## ---- results = 'hide', message = FALSE-------------------------------------------------
bform <- bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE)
fit_imp2 <- brm(bform, data = nhanes)

## ---------------------------------------------------------------------------------------
summary(fit_imp2)
conditional_effects(fit_imp2, "age:chl", resp = "bmi")

## ---------------------------------------------------------------------------------------
nhanes$se <- rexp(nrow(nhanes), 2)

## ---- results = 'hide', message = FALSE, eval = FALSE-----------------------------------
#  bform <- bf(bmi | mi() ~ age * mi(chl)) +
#    bf(chl | mi(se) ~ age) + set_rescor(FALSE)
#  fit_imp3 <- brm(bform, data = nhanes)
