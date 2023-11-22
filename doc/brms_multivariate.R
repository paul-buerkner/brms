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

## ----data-------------------------------------------------------------------------------
data("BTdata", package = "MCMCglmm")
head(BTdata)

## ----fit1, message=FALSE, warning=FALSE, results='hide'---------------------------------
bform1 <- 
  bf(mvbind(tarsus, back) ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam)) +
  set_rescor(TRUE)

fit1 <- brm(bform1, data = BTdata, chains = 2, cores = 2)

## ----summary1, warning=FALSE------------------------------------------------------------
fit1 <- add_criterion(fit1, "loo")
summary(fit1)

## ----pp_check1, message=FALSE-----------------------------------------------------------
pp_check(fit1, resp = "tarsus")
pp_check(fit1, resp = "back")

## ----R2_1-------------------------------------------------------------------------------
bayes_R2(fit1)

## ----fit2, message=FALSE, warning=FALSE, results='hide'---------------------------------
bf_tarsus <- bf(tarsus ~ sex + (1|p|fosternest) + (1|q|dam))
bf_back <- bf(back ~ hatchdate + (1|p|fosternest) + (1|q|dam))
fit2 <- brm(bf_tarsus + bf_back + set_rescor(TRUE), 
            data = BTdata, chains = 2, cores = 2)

## ----summary2, warning=FALSE------------------------------------------------------------
fit2 <- add_criterion(fit2, "loo")
summary(fit2)

## ----loo12------------------------------------------------------------------------------
loo(fit1, fit2)

## ----fit3, message=FALSE, warning=FALSE, results='hide'---------------------------------
bf_tarsus <- bf(tarsus ~ sex + (1|p|fosternest) + (1|q|dam)) +
  lf(sigma ~ 0 + sex) + skew_normal()
bf_back <- bf(back ~ s(hatchdate) + (1|p|fosternest) + (1|q|dam)) +
  gaussian()

fit3 <- brm(
  bf_tarsus + bf_back + set_rescor(FALSE),
  data = BTdata, chains = 2, cores = 2,
  control = list(adapt_delta = 0.95)
)

## ----summary3, warning=FALSE------------------------------------------------------------
fit3 <- add_criterion(fit3, "loo")
summary(fit3)

## ----me3--------------------------------------------------------------------------------
conditional_effects(fit3, "hatchdate", resp = "back")

