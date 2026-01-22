# Example: Ordered Beta Regression in brms
# Based on Kubinec (2023) - doi:10.1017/pan.2022.20

# install from my own working branch using remotes (if not installed already)

remotes::install_github("saudiwin/brms",ref="7e83794b3f56359b8d0799d9d0e5b41fdee9bbe2")

library(brms)

# Generate example data with 0s, 1s, and continuous values in (0,1)
set.seed(123)
n <- 200
x <- rnorm(n)

# True parameters
mu_true <- 0.5 * x
thres1_true <- -1   # threshold for boundary at 0 (parameter zoi in brms)
thres2_true <- 1.5  # threshold for boundary at 1 (parameter kappa = zoi + thres2_true in brms)

# Generate ordered beta data
pr_zero <- plogis(thres1_true - mu_true)
pr_one <- 1 - plogis(thres2_true - mu_true)
pr_cont <- plogis(thres2_true - mu_true) - plogis(thres1_true - mu_true)

y <- sapply(1:n, function(i) {
  component <- sample(1:3, 1, prob = c(pr_zero[i], pr_cont[i], pr_one[i]))
  if (component == 1) return(0)
  if (component == 3) return(1)
  rbeta(1, plogis(mu_true[i]) * 5, (1 - plogis(mu_true[i])) * 5)
})

dat <- data.frame(y = y, x = x)

# Check the data distribution
cat("Data summary:\n")
cat("  Zeros:", sum(dat$y == 0), "\n")
cat("  Ones:", sum(dat$y == 1), "\n")
cat("  Continuous (0,1):", sum(dat$y > 0 & dat$y < 1), "\n")

# Fit the ordered beta regression model
# The thresholds are estimated as ordered Intercepts (like ordinal models)
fit <- brm(
  y ~ x,
  data = dat,
  family = ordbeta(),
  chains = 2,
  cores = 2,
  iter = 2000
)

# View results
summary(fit)

# Posterior predictive check
# Note: this function doesn't work as it should due to discrete/continuous outcome
# makes fit look worse than it in fact is
pp_check(fit)

# Conditional effects
conditional_effects(fit)

# Get expected values
epred <- posterior_epred(fit)
cat("\nPosterior expected values dimensions:", dim(epred), "\n")

# Get predictions (including 0s and 1s)
pred <- posterior_predict(fit)
cat("Posterior predictions dimensions:", dim(pred), "\n")

# LOO cross-validation
loo_result <- loo(fit)
print(loo_result)

# -------------------------------------------------------------------
# Alternative: Model with phi (dispersion) varying
# -------------------------------------------------------------------

# fit_phi <- brm(
#   bf(y ~ x, phi ~ 1),
#   data = dat,
#   family = ordbeta(),
#   chains = 2,
#   cores = 2
# )
# summary(fit_phi)

# -------------------------------------------------------------------
# Just get the Stan code (without fitting)
# -------------------------------------------------------------------

cat("\n\nGenerated Stan code:\n")
cat("====================\n")
stancode(y ~ x, data = dat, family = ordbeta())
