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
#> Chain 1: Gradient evaluation took 0.000255 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.55 seconds.
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
#> Chain 1:  Elapsed Time: 4.181 seconds (Warm-up)
#> Chain 1:                2.965 seconds (Sampling)
#> Chain 1:                7.146 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000188 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.88 seconds.
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
#> Chain 2:  Elapsed Time: 4.141 seconds (Warm-up)
#> Chain 2:                2.962 seconds (Sampling)
#> Chain 2:                7.103 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000187 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.87 seconds.
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
#> Chain 3:  Elapsed Time: 4.245 seconds (Warm-up)
#> Chain 3:                2.96 seconds (Sampling)
#> Chain 3:                7.205 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.000179 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 1.79 seconds.
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
#> Chain 4:  Elapsed Time: 4.224 seconds (Warm-up)
#> Chain 4:                2.97 seconds (Sampling)
#> Chain 4:                7.194 seconds (Total)
#> Chain 4: 

# extract posterior samples of population-level effects
samples1 <- posterior_samples(fit, pars = "^b")
#> Warning: Method 'posterior_samples' is deprecated. Please see ?as_draws for recommended alternatives.
head(samples1)
#>   b_Intercept[1] b_Intercept[2] b_Intercept[3]    b_treat     b_period
#> 1      0.6342171       3.407759       4.924798 -1.2274213  0.238568168
#> 2      0.7047514       4.184251       5.218924 -0.9563724  0.269829086
#> 3      0.6502745       3.722971       5.410921 -0.5231338  0.268173971
#> 4      0.7485675       4.255780       5.103656 -1.3328422  0.451319378
#> 5      0.7370135       3.558930       5.094945 -0.9954435 -0.046348482
#> 6      0.7413192       3.965001       5.474634 -1.0588954 -0.008753989
#>       b_carry
#> 1  0.15806054
#> 2 -0.05294319
#> 3 -0.51239519
#> 4  0.10237665
#> 5 -0.36615368
#> 6 -0.25920380

# extract posterior samples of group-level standard deviations
samples2 <- posterior_samples(fit, pars = "^sd_")
#> Warning: Method 'posterior_samples' is deprecated. Please see ?as_draws for recommended alternatives.
head(samples2)
#>   sd_subject__Intercept
#> 1             0.9596501
#> 2             1.4156071
#> 3             1.2941843
#> 4             1.2327454
#> 5             1.3554111
#> 6             1.3621803
# }
```
