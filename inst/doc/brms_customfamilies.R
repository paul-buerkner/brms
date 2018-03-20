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
library(brms)
theme_set(theme_default())

## ----cbpp-------------------------------------------------------------------------------
data("cbpp", package = "lme4")
head(cbpp)

## ----fit1, results='hide'---------------------------------------------------------------
fit1 <- brm(incidence | trials(size) ~ period + (1|herd), 
            data = cbpp, family = binomial())

## ----fit1_summary-----------------------------------------------------------------------
summary(fit1)

## ----beta_binomial2---------------------------------------------------------------------
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "trials[n]"
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
stanvars <- stanvar(as.integer(cbpp$size), name = "trials")

## ----fit2, results='hide'---------------------------------------------------------------
fit2 <- brm(
  incidence ~ period + (1|herd), data = cbpp, 
  family = beta_binomial2, stan_funs = stan_funs,
  stanvars = stanvars
)

## ----summary_fit2-----------------------------------------------------------------------
summary(fit2)

## ---------------------------------------------------------------------------------------
expose_functions(fit2, vectorize = TRUE)

## ----log_lik----------------------------------------------------------------------------
log_lik_beta_binomial2 <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  N <- draws$data$trials[i]
  y <- draws$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, N)
}

## ----loo--------------------------------------------------------------------------------
LOO(fit1, fit2)

## ----predict----------------------------------------------------------------------------
predict_beta_binomial2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  N <- draws$data$trials[i]
  beta_binomial2_rng(mu, phi, N)
}

## ----pp_check---------------------------------------------------------------------------
pp_check(fit2)

## ----fitted-----------------------------------------------------------------------------
fitted_beta_binomial2 <- function(draws) {
  mu <- draws$dpars$mu
  trials <- draws$data$trials
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

## ----marginal_effects-------------------------------------------------------------------
marginal_effects(fit2, new_objects = list(trials = 1))

