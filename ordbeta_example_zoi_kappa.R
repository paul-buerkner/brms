# Example: Ordered Beta Regression in brms (with xi and kappa modeled)
# Based on Kubinec (2023) - doi:10.1017/pan.2022.20
# This version simulates covariate effects on BOTH xi and kappa thresholds
# and estimates them via brms distributional formulas.

# install from my own working branch using remotes (if not installed already)
remotes::install_github("saudiwin/brms", ref="4b7662bf51c4c4c5a39710fac3866f4efa634c3d")

library(brms)

# -------------------------------------------------------------------
# Simulate ordered beta data where xi and kappa depend on x
# -------------------------------------------------------------------
set.seed(123)
n <- 200
x <- rnorm(n)

# Mean model (mu) on the latent scale used by ordbeta
b_mu0_true <- 0.0
b_mux_true <- 0.5
mu_true <- b_mu0_true + b_mux_true * x

# xi model (lower threshold)
b_xi0_true <- -1.0
b_xix_true <- 0.8
xi_true <- b_xi0_true + b_xix_true * x

# kappa model (upper threshold), constructed to stay > xi
b_kappa0_true <- 1.0
b_kappa_true <- 0.3
kappa_true <- exp(b_kappa0_true + b_kappa_true * x) + xi_true # to stay positive

# Mixture probabilities for 0 / continuous / 1
# pr_zero = logistic(xi - mu)
# pr_one  = 1 - logistic(kappa - mu)
pr_zero <- plogis(xi_true - mu_true)
pr_one  <- 1 - plogis(kappa_true - mu_true)
pr_cont <- pmax(0, 1 - pr_zero - pr_one)

# Beta precision for the continuous part (kept constant for simplicity)
phi_true <- 5

y <- vapply(seq_len(n), function(i) {
  component <- sample(1:3, 1, prob = c(pr_zero[i], pr_cont[i], pr_one[i]))
  if (component == 1) return(0)
  if (component == 3) return(1)
  mu01 <- plogis(mu_true[i])
  rbeta(1, mu01 * phi_true, (1 - mu01) * phi_true)
}, numeric(1))

dat <- data.frame(y = y, x = x)

# Quick data check
cat("Data summary:\n")
cat("  Zeros:", sum(dat$y == 0), "\n")
cat("  Ones:", sum(dat$y == 1), "\n")
cat("  Continuous (0,1):", sum(dat$y > 0 & dat$y < 1), "\n")

# -------------------------------------------------------------------
# Fit ordered beta regression with distributional models for xi and kappa
# -------------------------------------------------------------------
# We estimate:
#   mu    ~ x
#   xi   ~ x
#   kappa ~ x
#
# NOTE: In the generative simulation above, kappa is constrained to exceed xi
# by construction. In estimation, brms will apply the ordbeta family's internal
# constraints/parameterization for ordered thresholds.

fit <- brm(
  bf(y ~ x, xi ~ x, kappa ~ x),
  data = dat,
  family = ordbeta(),
  chains = 2,
  cores = 2,
  iter = 2000
)

# View results
summary(fit)

# Posterior predictive check
pp_check(fit)

# Conditional effects for mu, xi, and kappa
conditional_effects(fit, dpar = "mu")
conditional_effects(fit, dpar = "xi")
conditional_effects(fit, dpar = "kappa")

# Expected values / predictions
epred <- posterior_epred(fit)
cat("\nPosterior expected values dimensions:", dim(epred), "\n")

pred <- posterior_predict(fit)
cat("Posterior predictions dimensions:", dim(pred), "\n")

# LOO cross-validation
loo_result <- loo(fit)
print(loo_result)

# -------------------------------------------------------------------
# Just get the Stan code (without fitting)
# -------------------------------------------------------------------
cat("\n\nGenerated Stan code:\n")
cat("====================\n")
stancode(bf(y ~ x, xi ~ x, kappa ~ x), data = dat, family = ordbeta())
