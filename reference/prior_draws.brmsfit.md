# Extract Prior Draws

Extract prior draws of specified parameters

## Usage

``` r
# S3 method for class 'brmsfit'
prior_draws(x, variable = NULL, pars = NULL, ...)

prior_draws(x, ...)

prior_samples(x, ...)
```

## Arguments

- x:

  An `R` object typically of class `brmsfit`.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

- pars:

  Deprecated alias of `variable`. For reasons of backwards
  compatibility, `pars` is interpreted as a vector of regular
  expressions by default unless `fixed = TRUE` is specified.

- ...:

  Arguments passed to individual methods (if applicable).

## Value

A `data.frame` containing the prior draws.

## Details

To make use of this function, the model must contain draws of prior
distributions. This can be ensured by setting `sample_prior = TRUE` in
function `brm`. Priors of certain parameters cannot be saved for
technical reasons. For instance, this is the case for the
population-level intercept, which is only computed after fitting the
model by default. If you want to treat the intercept as part of all the
other regression coefficients, so that sampling from its prior becomes
possible, use `... ~ 0 + Intercept + ...` in the formulas.

## Examples

``` r
# \dontrun{
fit <- brm(rating ~ treat + period + carry + (1|subject),
           data = inhaler, family = "cumulative",
           prior = set_prior("normal(0,2)", class = "b"),
           sample_prior = TRUE)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000365 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.65 seconds.
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
#> Chain 1:  Elapsed Time: 5.34 seconds (Warm-up)
#> Chain 1:                3.84 seconds (Sampling)
#> Chain 1:                9.18 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000235 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.35 seconds.
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
#> Chain 2:  Elapsed Time: 5.527 seconds (Warm-up)
#> Chain 2:                3.844 seconds (Sampling)
#> Chain 2:                9.371 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000237 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 2.37 seconds.
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
#> Chain 3:  Elapsed Time: 5.446 seconds (Warm-up)
#> Chain 3:                3.886 seconds (Sampling)
#> Chain 3:                9.332 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.000259 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 2.59 seconds.
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
#> Chain 4:  Elapsed Time: 5.526 seconds (Warm-up)
#> Chain 4:                3.833 seconds (Sampling)
#> Chain 4:                9.359 seconds (Total)
#> Chain 4: 

# extract all prior draws
draws1 <- prior_draws(fit)
head(draws1)
#>    Intercept           b  sd_subject
#> 1 -0.6914414  1.03966220  0.04305431
#> 2  0.4100858 -1.30219209  0.88032498
#> 3 -0.6499643 -4.22614007  4.23688415
#> 4  1.3236922 -2.28570696  3.99322064
#> 5  5.6541591 -0.09990052 14.30585129
#> 6 -8.0518432  3.43212761  3.95786077

# extract prior draws for the coefficient of 'treat'
draws2 <- prior_draws(fit, "b_treat")
head(draws2)
#>      b_treat
#> 1 -2.1192972
#> 2  1.5051579
#> 3 -0.5706313
#> 4 -4.8309096
#> 5  6.7716104
#> 6 -0.6459806
# }
```
