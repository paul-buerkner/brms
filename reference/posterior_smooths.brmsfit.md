# Posterior Predictions of Smooth Terms

Compute posterior predictions of smooth `s` and `t2` terms of models
fitted with brms.

## Usage

``` r
# S3 method for class 'brmsfit'
posterior_smooths(
  object,
  smooth,
  newdata = NULL,
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)

posterior_smooths(object, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- smooth:

  Name of a single smooth term for which predictions should be computed.

- newdata:

  An optional `data.frame` for which to evaluate predictions. If `NULL`
  (default), the original data of the model is used. Only those
  variables appearing in the chosen `smooth` term are required.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- dpar:

  Optional name of a predicted distributional parameter. If specified,
  expected predictions of this parameters are returned.

- nlpar:

  Optional name of a predicted non-linear parameter. If specified,
  expected predictions of this parameters are returned.

- ndraws:

  Positive integer indicating how many posterior draws should be used.
  If `NULL` (the default) all draws are used. Ignored if `draw_ids` is
  not `NULL`.

- draw_ids:

  An integer vector specifying the posterior draws to be used. If `NULL`
  (the default), all draws are used.

- ...:

  Currently ignored.

## Value

An S x N matrix, where S is the number of posterior draws and N is the
number of observations.

## Examples

``` r
# \dontrun{
set.seed(0)
dat <- mgcv::gamSim(1, n = 200, scale = 2)
#> Gu & Wahba 4 term additive model
fit <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.53 seconds.
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
#> Chain 1:  Elapsed Time: 5.796 seconds (Warm-up)
#> Chain 1:                4.633 seconds (Sampling)
#> Chain 1:                10.429 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.4 seconds.
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
#> Chain 2:  Elapsed Time: 5.771 seconds (Warm-up)
#> Chain 2:                4.819 seconds (Sampling)
#> Chain 2:                10.59 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4.8e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.48 seconds.
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
#> Chain 3:  Elapsed Time: 5.745 seconds (Warm-up)
#> Chain 3:                4.786 seconds (Sampling)
#> Chain 3:                10.531 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3.9e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 4:  Elapsed Time: 5.361 seconds (Warm-up)
#> Chain 4:                4.781 seconds (Sampling)
#> Chain 4:                10.142 seconds (Total)
#> Chain 4: 
#> Warning: There were 6 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit)
#> Warning: There were 6 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ s(x0) + s(x1) + s(x2) + s(x3) 
#>    Data: dat (Number of observations: 200) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Smoothing Spline Hyperparameters:
#>            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sds(sx0_1)     2.85      1.54     0.88     6.51 1.00     2217     1677
#> sds(sx1_1)     2.79      1.75     0.63     7.43 1.00     2301     2031
#> sds(sx2_1)    21.93      6.34    12.85    37.94 1.00     1451     2206
#> sds(sx3_1)     1.42      1.36     0.04     4.80 1.00     1725     1545
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     7.32      0.15     7.03     7.61 1.00     9272     2962
#> sx0_1        -1.40      5.67   -13.11     9.69 1.00     2894     2601
#> sx1_1        13.25      5.79     1.80    25.99 1.00     2492     2278
#> sx2_1        33.49     18.47    -3.80    69.18 1.00     2189     2304
#> sx3_1         0.57      3.90    -7.72     8.59 1.00     2181     1621
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.18      0.11     1.97     2.42 1.00     5664     2779
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

newdata <- data.frame(x2 = seq(0, 1, 10))
str(posterior_smooths(fit, smooth = "s(x2)", newdata = newdata))
#>  num [1:4000, 1] -4.66 -4.22 -5.15 -4.07 -4.57 ...
# }
```
