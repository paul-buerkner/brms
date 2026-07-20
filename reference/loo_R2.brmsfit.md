# Compute a LOO-adjusted R-squared for regression models

Compute a LOO-adjusted R-squared for regression models

## Usage

``` r
# S3 method for class 'brmsfit'
loo_R2(
  object,
  resp = NULL,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  seed = NULL,
  args_epred = list(),
  args_loglik = list(),
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead. Only used if `summary` is `TRUE`.

- probs:

  The percentiles to be computed by the `quantile` function. Only used
  if `summary` is `TRUE`.

- seed:

  Optional integer used to initialize the random number generator.

- args_epred:

  A named list of further arguments passed only to
  [`posterior_epred`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md).

- args_loglik:

  A named list of further arguments passed only to
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md).

- ...:

  Further arguments passed to both
  [`posterior_epred`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md)
  and
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md).

## Value

If `summary = TRUE`, an M x C matrix is returned (M = number of response
variables and c = `length(probs) + 2`) containing Bayesian bootstrap
based summary statistics of the LOO-adjusted R-squared values. If
`summary = FALSE`, the Bayesian bootstrap draws of the LOO-adjusted
R-squared values are returned in an S x M matrix (S is the number of
draws).

@details LOO-R2 uses LOO residuals and is defined as \\1-Var\_{loo-res}
/ Var_y\\, with \$\$ Var_y = V\_{n=1}^N y_n, and Var\_{loo-res} =
V\_{n=1}^N \hat{e}\_{loo,n}, \$\$ where
\\\hat{e}\_{loo,n}=y_n-\hat{y}\_{loo,n}\\. Bayesian bootstrap is used to
draw from the approximated uncertainty distribution as described by
Vehtari and Lampinen (2002).

## References

Vehtari and Lampinen (2002). Bayesian model assessment and comparison
using cross-validation predictive densities. Neural Computation,
14(10):2439-2468.

## Examples

``` r
# \dontrun{
fit <- brm(mpg ~ wt + cyl, data = mtcars)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
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
#> Chain 1:  Elapsed Time: 0.022 seconds (Warm-up)
#> Chain 1:                0.021 seconds (Sampling)
#> Chain 1:                0.043 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 2:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 2:                0.018 seconds (Sampling)
#> Chain 2:                0.039 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 3e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 3:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 3:                0.018 seconds (Sampling)
#> Chain 3:                0.039 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 4:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 4:                0.017 seconds (Sampling)
#> Chain 4:                0.038 seconds (Total)
#> Chain 4: 
summary(fit)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: mpg ~ wt + cyl 
#>    Data: mtcars (Number of observations: 32) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept    39.67      1.76    36.30    43.20 1.00     4766     2850
#> wt           -3.19      0.79    -4.74    -1.71 1.00     1909     2193
#> cyl          -1.51      0.43    -2.38    -0.64 1.00     1954     2082
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.66      0.36     2.07     3.51 1.00     2632     2341
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
loo_R2(fit)
#>     Estimate  Est.Error      Q2.5     Q97.5
#> R2 0.7854705 0.05390696 0.6605627 0.8764997

# compute R2 with new data
nd <- data.frame(mpg = c(10, 20, 30), wt = c(4, 3, 2), cyl = c(8, 6, 4))
loo_R2(fit, newdata = nd)
#>     Estimate    Est.Error      Q2.5     Q97.5
#> R2 0.8331886 3.595031e-05 0.8331055 0.8332681
# }
```
