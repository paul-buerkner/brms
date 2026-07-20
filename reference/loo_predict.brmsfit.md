# Compute Weighted Expectations Using LOO

These functions are wrappers around the
[`E_loo`](https://mc-stan.org/loo/reference/E_loo.html) function of the
loo package.

## Usage

``` r
# S3 method for class 'brmsfit'
loo_predict(
  object,
  type = c("mean", "var", "quantile"),
  probs = 0.5,
  psis_object = NULL,
  resp = NULL,
  ...
)

# S3 method for class 'brmsfit'
loo_epred(
  object,
  type = c("mean", "var", "quantile"),
  probs = 0.5,
  psis_object = NULL,
  resp = NULL,
  ...
)

loo_epred(object, ...)

# S3 method for class 'brmsfit'
loo_linpred(
  object,
  type = c("mean", "var", "quantile"),
  probs = 0.5,
  psis_object = NULL,
  resp = NULL,
  ...
)

# S3 method for class 'brmsfit'
loo_predictive_interval(object, prob = 0.9, psis_object = NULL, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- type:

  The statistic to be computed on the results. Can by either `"mean"`
  (default), `"var"`, or `"quantile"`.

- probs:

  A vector of quantiles to compute. Only used if `type = quantile`.

- psis_object:

  An optional object returned by
  [`psis`](https://mc-stan.org/loo/reference/psis.html). If
  `psis_object` is missing then
  [`psis`](https://mc-stan.org/loo/reference/psis.html) is executed
  internally, which may be time consuming for models fit to very large
  datasets.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- ...:

  Optional arguments passed to the underlying methods that is
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md),
  as well as
  [`posterior_predict`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md),
  [`posterior_epred`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md)
  or
  [`posterior_linpred`](https://paulbuerkner.com/brms/reference/posterior_linpred.brmsfit.md).

- prob:

  For `loo_predictive_interval`, a scalar in \\(0,1)\\ indicating the
  desired probability mass to include in the intervals. The default is
  `prob = 0.9` (\\90\\% intervals).

## Value

`loo_predict`, `loo_epred`, `loo_linpred`, and `loo_predictive_interval`
all return a matrix with one row per observation and one column per
summary statistic as specified by arguments `type` and `probs`. In
multivariate or categorical models a third dimension is added to
represent the response variables or categories, respectively.

`loo_predictive_interval(..., prob = p)` is equivalent to
`loo_predict(..., type = "quantile", probs = c(a, 1-a))` with
`a = (1 - p)/2`.

## Examples

``` r
# \dontrun{
## data from help("lm")
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
d <- data.frame(
  weight = c(ctl, trt),
  group = gl(2, 10, 20, labels = c("Ctl", "Trt"))
)
fit <- brm(weight ~ group, data = d)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
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
#> Chain 1:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 1:                0.011 seconds (Sampling)
#> Chain 1:                0.024 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 2:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 2:                0.012 seconds (Sampling)
#> Chain 2:                0.025 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 3:                0.013 seconds (Sampling)
#> Chain 3:                0.026 seconds (Total)
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
#> Chain 4:  Elapsed Time: 0.012 seconds (Warm-up)
#> Chain 4:                0.011 seconds (Sampling)
#> Chain 4:                0.023 seconds (Total)
#> Chain 4: 
loo_predictive_interval(fit, prob = 0.8)
#> Running PSIS to compute weights
#>            q10      q90
#>  [1,] 4.126218 6.126701
#>  [2,] 3.985866 5.946204
#>  [3,] 3.993827 6.084001
#>  [4,] 3.961473 5.841640
#>  [5,] 4.070708 6.100966
#>  [6,] 4.088472 6.125712
#>  [7,] 3.951910 6.049752
#>  [8,] 4.085542 6.063523
#>  [9,] 3.990548 6.029779
#> [10,] 3.998127 6.039939
#> [11,] 3.596160 5.691108
#> [12,] 3.705643 5.731761
#> [13,] 3.664478 5.717991
#> [14,] 3.859748 5.751465
#> [15,] 3.584900 5.397623
#> [16,] 3.764695 5.734718
#> [17,] 3.496421 5.378885
#> [18,] 3.578552 5.642172
#> [19,] 3.679650 5.707562
#> [20,] 3.648823 5.700667

## optionally log-weights can be pre-computed and reused
psis <- loo::psis(-log_lik(fit), cores = 2)
loo_predictive_interval(fit, prob = 0.8, psis_object = psis)
#>            q10      q90
#>  [1,] 4.175378 6.174271
#>  [2,] 3.970749 5.973336
#>  [3,] 3.981176 6.054476
#>  [4,] 3.982991 5.863508
#>  [5,] 4.079057 6.106985
#>  [6,] 4.038591 6.088208
#>  [7,] 3.984427 6.023762
#>  [8,] 4.060902 6.121451
#>  [9,] 3.989007 6.067537
#> [10,] 4.006101 6.026613
#> [11,] 3.590730 5.663438
#> [12,] 3.732675 5.755949
#> [13,] 3.591027 5.722428
#> [14,] 3.832287 5.738471
#> [15,] 3.618076 5.483889
#> [16,] 3.743359 5.730210
#> [17,] 3.618801 5.421882
#> [18,] 3.599374 5.640444
#> [19,] 3.687075 5.669565
#> [20,] 3.586922 5.659834
loo_predict(fit, type = "var", psis_object = psis)
#>             var
#>  [1,] 0.6480864
#>  [2,] 0.6580579
#>  [3,] 0.6609033
#>  [4,] 0.5635170
#>  [5,] 0.6464689
#>  [6,] 0.6690678
#>  [7,] 0.6626031
#>  [8,] 0.6602020
#>  [9,] 0.6748689
#> [10,] 0.6813778
#> [11,] 0.6519507
#> [12,] 0.6425287
#> [13,] 0.6835890
#> [14,] 0.6016842
#> [15,] 0.5486358
#> [16,] 0.6100506
#> [17,] 0.5641681
#> [18,] 0.6784335
#> [19,] 0.6646538
#> [20,] 0.6950113
loo_epred(fit, type = "var", psis_object = psis)
#>              var
#>  [1,] 0.06495195
#>  [2,] 0.06618208
#>  [3,] 0.06967061
#>  [4,] 0.05774600
#>  [5,] 0.06878128
#>  [6,] 0.06957301
#>  [7,] 0.06971152
#>  [8,] 0.06901383
#>  [9,] 0.06870589
#> [10,] 0.06982813
#> [11,] 0.06688398
#> [12,] 0.06495174
#> [13,] 0.06623014
#> [14,] 0.05870818
#> [15,] 0.05919239
#> [16,] 0.06180251
#> [17,] 0.05707503
#> [18,] 0.06677286
#> [19,] 0.06581242
#> [20,] 0.06688917
# }
```
