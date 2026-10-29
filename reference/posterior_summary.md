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
#> Chain 1:  Elapsed Time: 0.085 seconds (Warm-up)
#> Chain 1:                0.034 seconds (Sampling)
#> Chain 1:                0.119 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.081 seconds (Warm-up)
#> Chain 2:                0.029 seconds (Sampling)
#> Chain 2:                0.11 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.076 seconds (Warm-up)
#> Chain 3:                0.032 seconds (Sampling)
#> Chain 3:                0.108 seconds (Total)
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
#> Chain 4:  Elapsed Time: 0.083 seconds (Warm-up)
#> Chain 4:                0.036 seconds (Sampling)
#> Chain 4:                0.119 seconds (Total)
#> Chain 4: 
posterior_summary(fit)
#>                     Estimate   Est.Error        Q2.5       Q97.5
#> b_Intercept       20.1811232  86.3298543 -145.044694  192.147252
#> b_age              0.7906004   1.8747428   -3.011851    4.430065
#> b_sexfemale      170.2265425 101.3013191  -31.377733  365.896036
#> b_age:sexfemale   -2.5910860   2.1920688   -6.859041    1.774756
#> sigma            128.8263005  10.4167495  110.432145  151.924925
#> Intercept         96.5958410  14.4184857   68.497854  125.120103
#> lprior           -12.3419183   0.3919928  -13.201534  -11.652244
#> lp__            -485.0950476   1.5991424 -489.041916 -483.009295

# recommended workflow using posterior
library(posterior)
draws <- as_draws_array(fit)
summarise_draws(draws, default_summary_measures())
#> # A tibble: 8 × 7
#>   variable            mean   median      sd     mad       q5     q95
#>   <chr>              <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
#> 1 b_Intercept       20.2     19.5    86.3    85.7   -120.     164.  
#> 2 b_age              0.791    0.774   1.87    1.86    -2.34     3.77
#> 3 b_sexfemale      170.     170.    101.    100.       0.304  334.  
#> 4 b_age:sexfemale   -2.59    -2.60    2.19    2.18    -6.10     1.00
#> 5 sigma            129.     128.     10.4    10.2    113.     147.  
#> 6 Intercept         96.6     96.4    14.4    14.0     72.8    121.  
#> 7 lprior           -12.3    -12.3     0.392   0.386  -13.0    -11.8 
#> 8 lp__            -485.    -485.      1.60    1.40  -488.    -483.  
# }
```
