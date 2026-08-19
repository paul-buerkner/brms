# Summarize Posterior draws

Summarizes posterior draws based on point estimates (mean or median),
estimation errors (SD or MAD) and quantiles. This function mainly exists
to retain backwards compatibility. It will eventually be replaced by
functions of the posterior package (see examples below).

## Usage

``` r
posterior_summary(x, ...)

# Default S3 method
posterior_summary(x, probs = c(0.025, 0.975), robust = FALSE, ...)

# S3 method for class 'brmsfit'
posterior_summary(
  x,
  pars = NA,
  variable = NULL,
  probs = c(0.025, 0.975),
  robust = FALSE,
  ...
)
```

## Arguments

- x:

  An R object.

- ...:

  More arguments passed to or from other methods.

- probs:

  The percentiles to be computed by the
  [`quantile`](https://rdrr.io/r/stats/quantile.html) function.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead.

- pars:

  Deprecated alias of `variable`. For reasons of backwards
  compatibility, `pars` is interpreted as a vector of regular
  expressions by default unless `fixed = TRUE` is specified.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

## Value

A matrix where rows indicate variables and columns indicate the summary
estimates.

## See also

[`summarize_draws`](https://mc-stan.org/posterior/reference/draws_summary.html)

## Examples

``` r
# \dontrun{
fit <- brm(time ~ age * sex, data = kidney)
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
#> Chain 1:  Elapsed Time: 0.097 seconds (Warm-up)
#> Chain 1:                0.03 seconds (Sampling)
#> Chain 1:                0.127 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.101 seconds (Warm-up)
#> Chain 2:                0.045 seconds (Sampling)
#> Chain 2:                0.146 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.093 seconds (Warm-up)
#> Chain 3:                0.032 seconds (Sampling)
#> Chain 3:                0.125 seconds (Total)
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
#> Chain 4:  Elapsed Time: 0.1 seconds (Warm-up)
#> Chain 4:                0.037 seconds (Sampling)
#> Chain 4:                0.137 seconds (Total)
#> Chain 4: 
posterior_summary(fit)
#>                     Estimate   Est.Error        Q2.5       Q97.5
#> b_Intercept       18.5005937  87.7179795 -152.248782  184.178764
#> b_age              0.8126862   1.8874627   -2.756836    4.487141
#> b_sexfemale      170.2873881 103.3346862  -27.120354  367.237210
#> b_age:sexfemale   -2.5858524   2.2207673   -6.859541    1.639408
#> sigma            128.8130098  10.7099134  110.125291  151.409339
#> Intercept         96.0940203  14.7152594   67.391819  124.716570
#> lprior           -12.3308598   0.4007541  -13.184889  -11.617766
#> lp__            -485.1741313   1.6258448 -489.164938 -483.037219

# recommended workflow using posterior
library(posterior)
draws <- as_draws_array(fit)
summarise_draws(draws, default_summary_measures())
#> # A tibble: 8 × 7
#>   variable            mean   median      sd     mad      q5      q95
#>   <chr>              <dbl>    <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
#> 1 b_Intercept       18.5     20.7    87.7    88.3   -125.    160.   
#> 2 b_age              0.813    0.748   1.89    1.88    -2.22    3.95 
#> 3 b_sexfemale      170.     170.    103.    103.       1.66  338.   
#> 4 b_age:sexfemale   -2.59    -2.55    2.22    2.20    -6.26    0.923
#> 5 sigma            129.     128.     10.7    10.7    113.    148.   
#> 6 Intercept         96.1     96.2    14.7    14.4     71.2   120.   
#> 7 lprior           -12.3    -12.3     0.401   0.388  -13.0   -11.7  
#> 8 lp__            -485.    -485.      1.63    1.42  -488.   -483.   
# }
```
