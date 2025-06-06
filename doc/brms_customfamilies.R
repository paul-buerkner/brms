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

## ----cbpp-------------------------------------------------------------------------------
data("cbpp", package = "lme4")
head(cbpp)

## ----fit1, results='hide'---------------------------------------------------------------
fit1 <- brm(incidence | trials(size) ~ period + (1 | herd),
  data = cbpp, family = binomial()
)

## ----fit1_summary-----------------------------------------------------------------------
summary(fit1)

## ----beta_binomial2---------------------------------------------------------------------
beta_binomial2 <- custom_family(
  "beta_binomial2",
  dpars = c("mu", "phi"),
  links = c("logit", "log"),
  lb = c(0, 0), ub = c(1, NA),
  type = "int", vars = "vint1[n]"
)

## ----stan_funs--------------------------------------------------------------------------
stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"

## ----stanvars---------------------------------------------------------------------------
stanvars <- stanvar(scode = stan_funs, block = "functions")

## ----fit2, results='hide'---------------------------------------------------------------
fit2 <- brm(
  incidence | vint(size) ~ period + (1 | herd),
  data = cbpp,
  family = beta_binomial2, stanvars = stanvars
)

## ----summary_fit2-----------------------------------------------------------------------
summary(fit2)

## ---------------------------------------------------------------------------------------
expose_functions(fit2, vectorize = TRUE)

## ----log_lik----------------------------------------------------------------------------
log_lik_beta_binomial2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

## ----loo--------------------------------------------------------------------------------
loo(fit1, fit2)

## ----posterior_predict------------------------------------------------------------------
posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

## ----pp_check---------------------------------------------------------------------------
pp_check(fit2)

## ----posterior_epred--------------------------------------------------------------------
posterior_epred_beta_binomial2 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

## ----conditional_effects----------------------------------------------------------------
conditional_effects(fit2, conditions = data.frame(size = 1))
