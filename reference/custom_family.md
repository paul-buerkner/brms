# Custom Families in brms Models

Define custom families (i.e. response distribution) for use in brms
models. It allows users to benefit from the modeling flexibility of
brms, while applying their self-defined likelihood functions. All of the
post-processing methods for `brmsfit` objects can be made compatible
with custom families. See
[`vignette("brms_customfamilies")`](https://paulbuerkner.com/brms/articles/brms_customfamilies.md)
for more details. For a list of built-in families see
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).

## Usage

``` r
custom_family(
  name,
  dpars = "mu",
  links = "identity",
  type = c("real", "int"),
  lb = NA,
  ub = NA,
  vars = NULL,
  loop = TRUE,
  specials = NULL,
  threshold = "flexible",
  log_lik = NULL,
  posterior_predict = NULL,
  posterior_epred = NULL,
  predict = NULL,
  fitted = NULL,
  env = parent.frame()
)
```

## Arguments

- name:

  Name of the custom family.

- dpars:

  Names of the distributional parameters of the family. One parameter
  must be named `"mu"` and the main formula of the model will correspond
  to that parameter.

- links:

  Names of the link functions of the distributional parameters.

- type:

  Indicates if the response distribution is continuous (`"real"`) or
  discrete (`"int"`). This controls if the corresponding density
  function will be named with `<name>_lpdf` or `<name>_lpmf`.

- lb:

  Vector of lower bounds of the distributional parameters. Defaults to
  `NA` that is no lower bound.

- ub:

  Vector of upper bounds of the distributional parameters. Defaults to
  `NA` that is no upper bound.

- vars:

  Names of variables that are part of the likelihood function without
  being distributional parameters. That is, `vars` can be used to pass
  data to the likelihood. Such arguments will be added to the list of
  function arguments at the end, after the distributional parameters.
  See [`stanvar`](https://paulbuerkner.com/brms/reference/stanvar.md)
  for details about adding self-defined data to the generated Stan
  model. Addition arguments `vreal` and `vint` may be used for this
  purpose as well (see Examples below). See also
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  and
  [`addition-terms`](https://paulbuerkner.com/brms/reference/addition-terms.md)
  for more details.

- loop:

  Logical; Should the likelihood be evaluated via a loop (`TRUE`; the
  default) over observations in Stan? If `FALSE`, the Stan code will be
  written in a vectorized manner over observations if possible.

- specials:

  A character vector of special options to enable for this custom
  family. Currently for internal use only.

- threshold:

  Optional threshold type for custom ordinal families. Ignored for
  non-ordinal families.

- log_lik:

  Optional function to compute log-likelihood values of the model in R.
  This is only relevant if one wants to ensure compatibility with method
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md).

- posterior_predict:

  Optional function to compute posterior prediction of the model in R.
  This is only relevant if one wants to ensure compatibility with method
  [`posterior_predict`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md).

- posterior_epred:

  Optional function to compute expected values of the posterior
  predictive distribution of the model in R. This is only relevant if
  one wants to ensure compatibility with method
  [`posterior_epred`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md).

- predict:

  Deprecated alias of \`posterior_predict\`.

- fitted:

  Deprecated alias of \`posterior_epred\`.

- env:

  An [`environment`](https://rdrr.io/r/base/environment.html) in which
  certain post-processing functions related to the custom family can be
  found, if there were not directly passed to `custom_family`. This is
  only relevant if one wants to ensure compatibility with the methods
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md),
  [`posterior_predict`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md),
  or
  [`posterior_epred`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md).
  By default, `env` is the environment from which `custom_family` is
  called.

## Value

An object of class `customfamily` inheriting from class
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).

## Details

The corresponding probability density or mass `Stan` functions need to
have the same name as the custom family. That is if a family is called
`myfamily`, then the Stan functions should be called `myfamily_lpdf` or
`myfamily_lpmf` depending on whether it defines a continuous or discrete
distribution.

## See also

[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md),
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`stanvar`](https://paulbuerkner.com/brms/reference/stanvar.md)

## Examples

``` r
# \dontrun{
## demonstrate how to fit a beta-binomial model
## generate some fake data
phi <- 0.7
n <- 300
z <- rnorm(n, sd = 0.2)
ntrials <- sample(1:10, n, replace = TRUE)
eta <- 1 + z
mu <- exp(eta) / (1 + exp(eta))
a <- mu * phi
b <- (1 - mu) * phi
p <- rbeta(n, a, b)
y <- rbinom(n, ntrials, p)
dat <- data.frame(y, z, ntrials)

# define a custom family
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

# define the corresponding Stan density function
stan_density <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int N) {
    return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
  }
"
stanvars <- stanvar(scode = stan_density, block = "functions")

# fit the model
fit <- brm(y | vint(ntrials) ~ z, data = dat,
           family = beta_binomial2, stanvars = stanvars)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000138 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.38 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.736 seconds (Warm-up)
#> Chain 1:                0.665 seconds (Sampling)
#> Chain 1:                1.401 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000137 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.37 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.751 seconds (Warm-up)
#> Chain 2:                0.718 seconds (Sampling)
#> Chain 2:                1.469 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000129 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.29 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.73 seconds (Warm-up)
#> Chain 3:                0.673 seconds (Sampling)
#> Chain 3:                1.403 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.000139 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 1.39 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.726 seconds (Warm-up)
#> Chain 4:                0.659 seconds (Sampling)
#> Chain 4:                1.385 seconds (Total)
#> Chain 4: 
summary(fit)
#>  Family: beta_binomial2 
#>   Links: mu = logit 
#> Formula: y | vint(ntrials) ~ z 
#>    Data: dat (Number of observations: 300) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     1.01      0.11     0.80     1.22 1.00     3326     2971
#> z             1.36      0.53     0.35     2.40 1.00     3788     2926
#> 
#> Further Distributional Parameters:
#>     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> phi     0.80      0.11     0.61     1.04 1.00     3669     3103
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).


# define a *vectorized* custom family (no loop over observations)
# notice also that 'vint' no longer has an observation index
beta_binomial2_vec <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1", loop = FALSE
)

# define the corresponding Stan density function
stan_density_vec <- "
  real beta_binomial2_lpmf(array[] int y, vector mu, real phi, array[] int N) {
    return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
  }
"
stanvars_vec <- stanvar(scode = stan_density_vec, block = "functions")

# fit the model
fit_vec <- brm(y | vint(ntrials) ~ z, data = dat,
           family = beta_binomial2_vec,
           stanvars = stanvars_vec)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000118 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.18 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.527 seconds (Warm-up)
#> Chain 1:                0.456 seconds (Sampling)
#> Chain 1:                0.983 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.0001 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.528 seconds (Warm-up)
#> Chain 2:                0.55 seconds (Sampling)
#> Chain 2:                1.078 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000107 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.07 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.528 seconds (Warm-up)
#> Chain 3:                0.546 seconds (Sampling)
#> Chain 3:                1.074 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.000138 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 1.38 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.531 seconds (Warm-up)
#> Chain 4:                0.523 seconds (Sampling)
#> Chain 4:                1.054 seconds (Total)
#> Chain 4: 
summary(fit_vec)
#>  Family: beta_binomial2 
#>   Links: mu = logit 
#> Formula: y | vint(ntrials) ~ z 
#>    Data: dat (Number of observations: 300) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     1.01      0.11     0.80     1.21 1.00     3426     3160
#> z             1.38      0.53     0.37     2.43 1.00     3309     2470
#> 
#> Further Distributional Parameters:
#>     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> phi     0.80      0.11     0.61     1.05 1.00     3353     2920
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
# }
```
