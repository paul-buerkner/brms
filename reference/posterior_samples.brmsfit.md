# (Deprecated) Extract Posterior Samples

Extract posterior samples of specified parameters. The
`posterior_samples` method is deprecated. We recommend using the more
modern and consistent
[`as_draws_*`](https://paulbuerkner.com/brms/reference/draws-brms.md)
extractor functions of the posterior package instead.

## Usage

``` r
# S3 method for class 'brmsfit'
posterior_samples(
  x,
  pars = NA,
  fixed = FALSE,
  add_chain = FALSE,
  subset = NULL,
  as.matrix = FALSE,
  as.array = FALSE,
  ...
)

posterior_samples(x, pars = NA, ...)
```

## Arguments

- x:

  An `R` object typically of class `brmsfit`

- pars:

  Names of parameters for which posterior samples should be returned, as
  given by a character vector or regular expressions. By default, all
  posterior samples of all parameters are extracted.

- fixed:

  Indicates whether parameter names should be matched exactly (`TRUE`)
  or treated as regular expressions (`FALSE`). Default is `FALSE`.

- add_chain:

  A flag indicating if the returned `data.frame` should contain two
  additional columns. The `chain` column indicates the chain in which
  each sample was generated, the `iter` column indicates the iteration
  number within each chain.

- subset:

  A numeric vector indicating the rows (i.e., posterior samples) to be
  returned. If `NULL` (the default), all posterior samples are returned.

- as.matrix:

  Should the output be a `matrix` instead of a `data.frame`? Defaults to
  `FALSE`.

- as.array:

  Should the output be an `array` instead of a `data.frame`? Defaults to
  `FALSE`.

- ...:

  Arguments passed to individual methods (if applicable).

## Value

A data.frame (matrix or array) containing the posterior samples.

## See also

[`as_draws`](https://paulbuerkner.com/brms/reference/draws-brms.md),
[`as.data.frame`](https://paulbuerkner.com/brms/reference/as.data.frame.brmsfit.md)

## Examples

``` r
# \dontrun{
fit <- brm(rating ~ treat + period + carry + (1|subject),
           data = inhaler, family = "cumulative")
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000256 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.56 seconds.
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
#> Chain 1:  Elapsed Time: 5.355 seconds (Warm-up)
#> Chain 1:                3.845 seconds (Sampling)
#> Chain 1:                9.2 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000238 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.38 seconds.
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
#> Chain 2:  Elapsed Time: 5.391 seconds (Warm-up)
#> Chain 2:                3.809 seconds (Sampling)
#> Chain 2:                9.2 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000238 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 2.38 seconds.
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
#> Chain 3:  Elapsed Time: 5.482 seconds (Warm-up)
#> Chain 3:                3.822 seconds (Sampling)
#> Chain 3:                9.304 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.000231 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 2.31 seconds.
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
#> Chain 4:  Elapsed Time: 5.756 seconds (Warm-up)
#> Chain 4:                3.813 seconds (Sampling)
#> Chain 4:                9.569 seconds (Total)
#> Chain 4: 

# extract posterior samples of population-level effects
samples1 <- posterior_samples(fit, pars = "^b")
#> Warning: Method 'posterior_samples' is deprecated. Please see ?as_draws for recommended alternatives.
head(samples1)
#>   b_Intercept[1] b_Intercept[2] b_Intercept[3]    b_treat   b_period
#> 1      0.6597909       4.049887       5.129223 -1.4129508 0.17261708
#> 2      0.7894877       3.763109       5.339839 -0.9981397 0.38029214
#> 3      0.7792166       4.146578       5.317411 -0.9564157 0.50993634
#> 4      0.6732261       3.988528       5.475101 -0.4825884 0.07534019
#> 5      0.7056669       4.340685       5.499632 -1.5816312 0.44556879
#> 6      0.4940245       3.781961       5.445265 -0.9741098 0.10051181
#>       b_carry
#> 1 -0.12440360
#> 2 -0.49657015
#> 3 -0.28599545
#> 4 -0.57014270
#> 5 -0.03337895
#> 6 -0.50983136

# extract posterior samples of group-level standard deviations
samples2 <- posterior_samples(fit, pars = "^sd_")
#> Warning: Method 'posterior_samples' is deprecated. Please see ?as_draws for recommended alternatives.
head(samples2)
#>   sd_subject__Intercept
#> 1              1.419671
#> 2              1.111153
#> 3              1.796667
#> 4              1.470518
#> 5              1.738034
#> 6              1.273054
# }
```
