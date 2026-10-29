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
#> Chain 1: Gradient evaluation took 5e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
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
#> Chain 1:  Elapsed Time: 0.011 seconds (Warm-up)
#> Chain 1:                0.011 seconds (Sampling)
#> Chain 1:                0.022 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.011 seconds (Warm-up)
#> Chain 2:                0.01 seconds (Sampling)
#> Chain 2:                0.021 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
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
#> Chain 3:  Elapsed Time: 0.01 seconds (Warm-up)
#> Chain 3:                0.011 seconds (Sampling)
#> Chain 3:                0.021 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
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
#> Chain 4:                0.01 seconds (Sampling)
#> Chain 4:                0.022 seconds (Total)
#> Chain 4: 
loo_predictive_interval(fit, prob = 0.8)
#> Running PSIS to compute weights
#>            q10      q90
#>  [1,] 4.104461 6.106387
#>  [2,] 3.946844 5.942279
#>  [3,] 3.975836 6.010287
#>  [4,] 3.936007 5.859922
#>  [5,] 4.104679 6.085987
#>  [6,] 4.082916 6.108275
#>  [7,] 3.965488 6.054453
#>  [8,] 4.056460 6.076927
#>  [9,] 3.985012 6.056993
#> [10,] 3.988782 6.066735
#> [11,] 3.626811 5.694474
#> [12,] 3.691668 5.733745
#> [13,] 3.655181 5.692952
#> [14,] 3.856260 5.747909
#> [15,] 3.611965 5.522781
#> [16,] 3.794947 5.688854
#> [17,] 3.594747 5.418993
#> [18,] 3.585373 5.706027
#> [19,] 3.663547 5.736795
#> [20,] 3.617765 5.715515

## optionally log-weights can be pre-computed and reused
psis <- loo::psis(-log_lik(fit), cores = 2)
loo_predictive_interval(fit, prob = 0.8, psis_object = psis)
#>            q10      q90
#>  [1,] 4.106413 6.072142
#>  [2,] 3.956705 5.988746
#>  [3,] 4.012618 6.042415
#>  [4,] 3.925902 5.843192
#>  [5,] 4.106890 6.122774
#>  [6,] 4.065625 6.124774
#>  [7,] 3.968849 5.993979
#>  [8,] 4.044554 6.095990
#>  [9,] 3.958221 6.024337
#> [10,] 4.042580 6.063134
#> [11,] 3.608261 5.666594
#> [12,] 3.694262 5.762264
#> [13,] 3.654366 5.744219
#> [14,] 3.814131 5.753113
#> [15,] 3.628184 5.476126
#> [16,] 3.773774 5.726626
#> [17,] 3.616112 5.404928
#> [18,] 3.601986 5.701068
#> [19,] 3.683811 5.692551
#> [20,] 3.622785 5.692261
loo_predict(fit, type = "var", psis_object = psis)
#>             var
#>  [1,] 0.6124240
#>  [2,] 0.6479519
#>  [3,] 0.6599563
#>  [4,] 0.5830894
#>  [5,] 0.6373206
#>  [6,] 0.6332037
#>  [7,] 0.6730528
#>  [8,] 0.6409808
#>  [9,] 0.6613779
#> [10,] 0.7089992
#> [11,] 0.6953620
#> [12,] 0.6129567
#> [13,] 0.6869191
#> [14,] 0.5858355
#> [15,] 0.5375803
#> [16,] 0.6318167
#> [17,] 0.5301938
#> [18,] 0.6679234
#> [19,] 0.6906654
#> [20,] 0.6373395
loo_epred(fit, type = "var", psis_object = psis)
#>              var
#>  [1,] 0.05941590
#>  [2,] 0.06079187
#>  [3,] 0.06297466
#>  [4,] 0.05487445
#>  [5,] 0.06222280
#>  [6,] 0.06274935
#>  [7,] 0.06302882
#>  [8,] 0.06243327
#>  [9,] 0.06238772
#> [10,] 0.06313856
#> [11,] 0.06661047
#> [12,] 0.06479759
#> [13,] 0.06630318
#> [14,] 0.05648581
#> [15,] 0.05430182
#> [16,] 0.06063892
#> [17,] 0.05064226
#> [18,] 0.06637826
#> [19,] 0.06587703
#> [20,] 0.06679601
# }
```
