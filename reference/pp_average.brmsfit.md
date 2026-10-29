# Posterior predictive draws averaged across models

Compute posterior predictive draws averaged across models. Weighting can
be done in various ways, for instance using Akaike weights based on
information criteria or marginal likelihoods.

## Usage

``` r
# S3 method for class 'brmsfit'
pp_average(
  x,
  ...,
  weights = "stacking",
  method = "posterior_predict",
  ndraws = NULL,
  nsamples = NULL,
  summary = TRUE,
  probs = c(0.025, 0.975),
  robust = FALSE,
  model_names = NULL,
  control = list(),
  seed = NULL
)

pp_average(x, ...)
```

## Arguments

- x:

  A `brmsfit` object.

- ...:

  More `brmsfit` objects or further arguments passed to the underlying
  post-processing functions. In particular, see
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  for further supported arguments.

- weights:

  Name of the criterion to compute weights from. Should be one of
  `"loo"`, `"waic"`, `"kfold"`, `"stacking"` (current default), `"bma"`,
  or `"pseudobma"`. For the former three options, Akaike weights will be
  computed based on the information criterion values returned by the
  respective methods. For `"stacking"` and `"pseudobma"`, method
  [`loo_model_weights`](https://paulbuerkner.com/brms/reference/loo_model_weights.brmsfit.md)
  will be used to obtain weights. For `"bma"`, method
  [`post_prob`](https://paulbuerkner.com/brms/reference/post_prob.brmsfit.md)
  will be used to compute Bayesian model averaging weights based on log
  marginal likelihood values (make sure to specify reasonable priors in
  this case). For some methods, `weights` may also be a numeric vector
  of pre-specified weights.

- method:

  Method used to obtain predictions to average over. Should be one of
  `"posterior_predict"` (default), `"posterior_epred"`,
  `"posterior_linpred"` or `"predictive_error"`.

- ndraws:

  Total number of posterior draws to use.

- nsamples:

  Deprecated alias of `ndraws`.

- summary:

  Should summary statistics (i.e. means, sds, and 95% intervals) be
  returned instead of the raw values? Default is `TRUE`.

- probs:

  The percentiles to be computed by the `quantile` function. Only used
  if `summary` is `TRUE`.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead. Only used if `summary` is `TRUE`.

- model_names:

  If `NULL` (the default) will use model names derived from deparsing
  the call. Otherwise will use the passed values as model names.

- control:

  Optional `list` of further arguments passed to the function specified
  in `weights`.

- seed:

  A single numeric value passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) to make results
  reproducible.

## Value

Same as the output of the method specified in argument `method`.

## Details

Weights are computed with the
[`model_weights`](https://paulbuerkner.com/brms/reference/model_weights.brmsfit.md)
method.

## See also

[`model_weights`](https://paulbuerkner.com/brms/reference/model_weights.brmsfit.md),
[`posterior_average`](https://paulbuerkner.com/brms/reference/posterior_average.brmsfit.md)

## Examples

``` r
# \dontrun{
# model with 'treat' as predictor
fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
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
#> Chain 1:  Elapsed Time: 0.023 seconds (Warm-up)
#> Chain 1:                0.021 seconds (Sampling)
#> Chain 1:                0.044 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.022 seconds (Warm-up)
#> Chain 2:                0.021 seconds (Sampling)
#> Chain 2:                0.043 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.023 seconds (Warm-up)
#> Chain 3:                0.02 seconds (Sampling)
#> Chain 3:                0.043 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 4:  Elapsed Time: 0.022 seconds (Warm-up)
#> Chain 4:                0.021 seconds (Sampling)
#> Chain 4:                0.043 seconds (Total)
#> Chain 4: 
summary(fit1)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: rating ~ treat + period + carry 
#>    Data: inhaler (Number of observations: 572) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     1.43      0.03     1.38     1.48 1.00     4589     2877
#> treat        -0.25      0.07    -0.39    -0.11 1.00     3616     2925
#> period        0.05      0.05    -0.05     0.15 1.00     4587     2860
#> carry        -0.04      0.05    -0.15     0.06 1.00     3538     3239
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     0.61      0.02     0.57     0.65 1.00     5194     2878
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

# model without 'treat' as predictor
fit2 <- brm(rating ~ period + carry, data = inhaler)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
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
#> Chain 1:  Elapsed Time: 0.018 seconds (Warm-up)
#> Chain 1:                0.017 seconds (Sampling)
#> Chain 1:                0.035 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.017 seconds (Warm-up)
#> Chain 2:                0.014 seconds (Sampling)
#> Chain 2:                0.031 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 3:  Elapsed Time: 0.018 seconds (Warm-up)
#> Chain 3:                0.015 seconds (Sampling)
#> Chain 3:                0.033 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 4:  Elapsed Time: 0.018 seconds (Warm-up)
#> Chain 4:                0.016 seconds (Sampling)
#> Chain 4:                0.034 seconds (Total)
#> Chain 4: 
summary(fit2)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: rating ~ period + carry 
#>    Data: inhaler (Number of observations: 572) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     1.43      0.02     1.38     1.48 1.00     4023     2800
#> period        0.05      0.05    -0.05     0.15 1.00     4763     3007
#> carry        -0.17      0.04    -0.24    -0.10 1.00     4795     3321
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     0.61      0.02     0.58     0.65 1.00     4074     2465
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

# compute model-averaged predicted values
(df <- unique(inhaler[, c("treat", "period", "carry")]))
#>     treat period carry
#> 1     0.5    0.5     0
#> 60   -0.5    0.5     0
#> 287  -0.5   -0.5    -1
#> 346   0.5   -0.5     1
pp_average(fit1, fit2, newdata = df)
#>      Estimate Est.Error       Q2.5    Q97.5
#> [1,] 1.355118 0.6052531 0.17856944 2.536315
#> [2,] 1.575359 0.5989236 0.42295301 2.744175
#> [3,] 1.596402 0.6044013 0.42720091 2.770042
#> [4,] 1.240386 0.6170943 0.01388588 2.432680
#> attr(,"weights")
#>      fit1      fit2 
#> 0.8516506 0.1483494 
#> attr(,"ndraws")
#> fit1 fit2 
#> 3407  593 

# compute model-averaged fitted values
pp_average(fit1, fit2, method = "fitted", newdata = df)
#>      Estimate  Est.Error     Q2.5    Q97.5
#> [1,] 1.349206 0.06765766 1.229548 1.491723
#> [2,] 1.564406 0.06526306 1.427076 1.677502
#> [3,] 1.577877 0.05157725 1.475464 1.677291
#> [4,] 1.235182 0.05214174 1.132670 1.335806
#> attr(,"weights")
#>      fit1      fit2 
#> 0.8516506 0.1483494 
#> attr(,"ndraws")
#> fit1 fit2 
#> 3407  593 
# }
```
