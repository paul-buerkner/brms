## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
#  b <- c(2, 0.75)
#  x <- rnorm(100)
#  y <- rnorm(100, mean = b[1] * exp(b[2] * x))
#  dat1 <- data.frame(x, y)

## ------------------------------------------------------------------------
#  library(brms)
#  prior1 <- c(prior(normal(1, 2), nlpar = "b1"),
#              prior(normal(0, 2), nlpar = "b2"))
#  fit1 <- brm(y ~ b1 * exp(b2 * x), nonlinear = b1 + b2 ~ 1,
#              data = dat1, prior = prior1)

## ------------------------------------------------------------------------
#  summary(fit1)
#  plot(fit1)
#  plot(marginal_effects(fit1), points = TRUE)

## ------------------------------------------------------------------------
#  fit2 <- brm(y ~ x, data = dat1)
#  summary(fit2)

## ------------------------------------------------------------------------
#  pp_check(fit1)
#  pp_check(fit2)

## ------------------------------------------------------------------------
#  LOO(fit1, fit2)

## ------------------------------------------------------------------------
#  url <- "https://raw.githubusercontent.com/mages/diesunddas/master/Data/ClarkTriangle.csv"
#  loss <- read.csv(url)
#  head(loss)

## ------------------------------------------------------------------------
#  fit_loss <- brm(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
#                  nonlinear = list(ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1),
#                  data = loss, family = gaussian(),
#                  prior = c(prior(normal(5000, 1000), nlpar = "ult"),
#                            prior(normal(1, 2), nlpar = "omega"),
#                            prior(normal(45, 10), nlpar = "theta")),
#                  control = list(adapt_delta = 0.9))

## ------------------------------------------------------------------------
#  summary(fit_loss)
#  plot(fit_loss)
#  marginal_effects(fit_loss)

## ------------------------------------------------------------------------
#  conditions <- data.frame(AY = unique(loss$AY))
#  rownames(conditions) <- unique(loss$AY)
#  plot(marginal_effects(fit_loss, conditions = conditions,
#                        re_formula = NULL, method = "predict"),
#       ncol = 5, points = TRUE)

## ------------------------------------------------------------------------
#  inv_logit <- function(x) 1 / (1 + exp(-x))
#  ability <- rnorm(300)
#  p <- 0.33 + 0.67 * inv_logit(ability)
#  answer <- ifelse(runif(300, 0, 1) < p, 1, 0)
#  dat_ir <- data.frame(ability, answer)

## ------------------------------------------------------------------------
#  fit_ir1 <- brm(answer ~ ability, data = dat_ir, family = bernoulli())

## ------------------------------------------------------------------------
#  summary(fit_ir1)
#  plot(marginal_effects(fit_ir1), points = TRUE)

## ------------------------------------------------------------------------
#  fit_ir2 <- brm(answer ~ 0.33 + 0.67 * inv_logit(eta),
#                 nonlinear = eta ~ ability,
#                 data = dat_ir, family = bernoulli("identity"),
#                 prior = prior(normal(0, 5), nlpar = "eta"))

## ------------------------------------------------------------------------
#  summary(fit_ir2)
#  plot(marginal_effects(fit_ir2), points = TRUE)

## ------------------------------------------------------------------------
#  LOO(fit_ir1, fit_ir2)

## ------------------------------------------------------------------------
#  fit_ir3 <- brm(answer ~ guess + (1 - guess) * inv_logit(eta),
#                 nonlinear = list(eta ~ 0 + ability, guess ~ 1),
#                 data = dat_ir, family = bernoulli("identity"),
#                 prior = c(prior(normal(0, 5), nlpar = "eta"),
#                           prior(beta(1, 1), nlpar = "guess", lb = 0, ub = 1)))

## ------------------------------------------------------------------------
#  summary(fit_ir3)
#  plot(fit_ir3)
#  plot(marginal_effects(fit_ir3), points = TRUE)

