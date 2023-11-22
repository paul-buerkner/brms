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
b <- c(2, 0.75)
x <- rnorm(100)
y <- rnorm(100, mean = b[1] * exp(b[2] * x))
dat1 <- data.frame(x, y)

## ---- results='hide'--------------------------------------------------------------------
prior1 <- prior(normal(1, 2), nlpar = "b1") +
  prior(normal(0, 2), nlpar = "b2")
fit1 <- brm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
            data = dat1, prior = prior1)

## ---------------------------------------------------------------------------------------
summary(fit1)
plot(fit1)
plot(conditional_effects(fit1), points = TRUE)

## ---- results='hide'--------------------------------------------------------------------
fit2 <- brm(y ~ x, data = dat1)

## ---------------------------------------------------------------------------------------
summary(fit2)

## ---------------------------------------------------------------------------------------
pp_check(fit1)
pp_check(fit2)

## ---------------------------------------------------------------------------------------
loo(fit1, fit2)

## ---------------------------------------------------------------------------------------
data(loss)
head(loss)

## ---- results='hide'--------------------------------------------------------------------
fit_loss <- brm(
  bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
     ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
     nl = TRUE),
  data = loss, family = gaussian(),
  prior = c(
    prior(normal(5000, 1000), nlpar = "ult"),
    prior(normal(1, 2), nlpar = "omega"),
    prior(normal(45, 10), nlpar = "theta")
  ),
  control = list(adapt_delta = 0.9)
)

## ---------------------------------------------------------------------------------------
summary(fit_loss)
plot(fit_loss, N = 3, ask = FALSE)
conditional_effects(fit_loss)

## ---------------------------------------------------------------------------------------
conditions <- data.frame(AY = unique(loss$AY))
rownames(conditions) <- unique(loss$AY)
me_loss <- conditional_effects(
  fit_loss, conditions = conditions,
  re_formula = NULL, method = "predict"
)
plot(me_loss, ncol = 5, points = TRUE)

## ---------------------------------------------------------------------------------------
inv_logit <- function(x) 1 / (1 + exp(-x))
ability <- rnorm(300)
p <- 0.33 + 0.67 * inv_logit(ability)
answer <- ifelse(runif(300, 0, 1) < p, 1, 0)
dat_ir <- data.frame(ability, answer)

## ---- results='hide'--------------------------------------------------------------------
fit_ir1 <- brm(answer ~ ability, data = dat_ir, family = bernoulli())

## ---------------------------------------------------------------------------------------
summary(fit_ir1)
plot(conditional_effects(fit_ir1), points = TRUE)

## ---- results='hide'--------------------------------------------------------------------
fit_ir2 <- brm(
  bf(answer ~ 0.33 + 0.67 * inv_logit(eta),
     eta ~ ability, nl = TRUE),
  data = dat_ir, family = bernoulli("identity"),
  prior = prior(normal(0, 5), nlpar = "eta")
)

## ---------------------------------------------------------------------------------------
summary(fit_ir2)
plot(conditional_effects(fit_ir2), points = TRUE)

## ---------------------------------------------------------------------------------------
loo(fit_ir1, fit_ir2)

## ---- results='hide'--------------------------------------------------------------------
fit_ir3 <- brm(
  bf(answer ~ guess + (1 - guess) * inv_logit(eta),
    eta ~ 0 + ability, guess ~ 1, nl = TRUE),
  data = dat_ir, family = bernoulli("identity"),
  prior = c(
    prior(normal(0, 5), nlpar = "eta"),
    prior(beta(1, 1), nlpar = "guess", lb = 0, ub = 1)
  )
)

## ---------------------------------------------------------------------------------------
summary(fit_ir3)
plot(fit_ir3)
plot(conditional_effects(fit_ir3), points = TRUE)

