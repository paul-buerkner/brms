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
#> Chain 1: Gradient evaluation took 9e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
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
#> Chain 1:  Elapsed Time: 0.012 seconds (Warm-up)
#> Chain 1:                0.011 seconds (Sampling)
#> Chain 1:                0.023 seconds (Total)
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
#> Chain 2:                0.011 seconds (Sampling)
#> Chain 2:                0.022 seconds (Total)
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
#> Chain 3:                0.012 seconds (Sampling)
#> Chain 3:                0.025 seconds (Total)
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
#> Chain 4:                0.01 seconds (Sampling)
#> Chain 4:                0.022 seconds (Total)
#> Chain 4: 
loo_predictive_interval(fit, prob = 0.8)
#> Running PSIS to compute weights
#>            q10      q90
#>  [1,] 4.086658 6.099511
#>  [2,] 3.952639 5.982989
#>  [3,] 3.988170 6.003492
#>  [4,] 3.989440 5.835315
#>  [5,] 4.086716 6.089946
#>  [6,] 4.056988 6.106231
#>  [7,] 3.950568 6.050696
#>  [8,] 4.090262 6.057129
#>  [9,] 3.973462 6.073383
#> [10,] 3.989869 6.054655
#> [11,] 3.641428 5.638272
#> [12,] 3.721266 5.721369
#> [13,] 3.643473 5.713548
#> [14,] 3.866520 5.761594
#> [15,] 3.591131 5.482531
#> [16,] 3.790882 5.747139
#> [17,] 3.592316 5.409894
#> [18,] 3.599975 5.728514
#> [19,] 3.686796 5.709781
#> [20,] 3.610006 5.758691

## optionally log-weights can be pre-computed and reused
psis <- loo::psis(-log_lik(fit), cores = 2)
loo_predictive_interval(fit, prob = 0.8, psis_object = psis)
#>            q10      q90
#>  [1,] 4.091342 6.130623
#>  [2,] 3.917652 6.000012
#>  [3,] 4.010131 6.035042
#>  [4,] 3.965405 5.837825
#>  [5,] 4.118871 6.110691
#>  [6,] 4.047368 6.102830
#>  [7,] 3.991166 5.986654
#>  [8,] 4.050299 6.097324
#>  [9,] 3.945654 6.006360
#> [10,] 4.017182 6.073838
#> [11,] 3.621182 5.671400
#> [12,] 3.738226 5.752884
#> [13,] 3.686083 5.765674
#> [14,] 3.828719 5.708753
#> [15,] 3.643854 5.455236
#> [16,] 3.782444 5.718142
#> [17,] 3.494880 5.364733
#> [18,] 3.653488 5.706662
#> [19,] 3.682481 5.738590
#> [20,] 3.640418 5.706361
loo_predict(fit, type = "var", psis_object = psis)
#>             var
#>  [1,] 0.6341889
#>  [2,] 0.6459611
#>  [3,] 0.6477019
#>  [4,] 0.5838173
#>  [5,] 0.6426555
#>  [6,] 0.6435834
#>  [7,] 0.6558330
#>  [8,] 0.6610992
#>  [9,] 0.6480890
#> [10,] 0.6780068
#> [11,] 0.6663567
#> [12,] 0.6302480
#> [13,] 0.6854037
#> [14,] 0.5926749
#> [15,] 0.5611214
#> [16,] 0.6322851
#> [17,] 0.5149118
#> [18,] 0.6851024
#> [19,] 0.6877240
#> [20,] 0.6617218
loo_epred(fit, type = "var", psis_object = psis)
#>              var
#>  [1,] 0.05823207
#>  [2,] 0.06063710
#>  [3,] 0.06278810
#>  [4,] 0.05323118
#>  [5,] 0.06129478
#>  [6,] 0.06195281
#>  [7,] 0.06280742
#>  [8,] 0.06149881
#>  [9,] 0.06226685
#> [10,] 0.06286526
#> [11,] 0.06773635
#> [12,] 0.06530733
#> [13,] 0.06716529
#> [14,] 0.05653938
#> [15,] 0.05727392
#> [16,] 0.06101407
#> [17,] 0.05522012
#> [18,] 0.06744473
#> [19,] 0.06657905
#> [20,] 0.06789751
# }
```
