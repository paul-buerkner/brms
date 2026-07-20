# Extract Diagnostic Quantities of brms Models

Extract quantities that can be used to diagnose sampling behavior of the
algorithms applied by Stan at the back-end of brms.

## Usage

``` r
# S3 method for class 'brmsfit'
log_posterior(object, ...)

# S3 method for class 'brmsfit'
nuts_params(object, pars = NULL, ...)

# S3 method for class 'brmsfit'
rhat(x, pars = NULL, ...)

# S3 method for class 'brmsfit'
neff_ratio(object, pars = NULL, ...)
```

## Arguments

- object, x:

  A `brmsfit` object.

- ...:

  Arguments passed to individual methods.

- pars:

  An optional character vector of parameter names. For `nuts_params`
  these will be NUTS sampler parameter names rather than model
  parameters. If pars is omitted all parameters are included.

## Value

The exact form of the output depends on the method.

## Details

For more details see
[`bayesplot-extractors`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html).

## Examples

``` r
# \dontrun{
fit <- brm(time ~ age * sex, data = kidney)
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
#> Chain 1:  Elapsed Time: 0.098 seconds (Warm-up)
#> Chain 1:                0.033 seconds (Sampling)
#> Chain 1:                0.131 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.116 seconds (Warm-up)
#> Chain 2:                0.042 seconds (Sampling)
#> Chain 2:                0.158 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.104 seconds (Warm-up)
#> Chain 3:                0.034 seconds (Sampling)
#> Chain 3:                0.138 seconds (Total)
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
#> Chain 4:  Elapsed Time: 0.108 seconds (Warm-up)
#> Chain 4:                0.034 seconds (Sampling)
#> Chain 4:                0.142 seconds (Total)
#> Chain 4: 

lp <- log_posterior(fit)
head(lp)
#>   Chain Iteration     Value
#> 1     1         1 -483.6940
#> 2     1         2 -483.4945
#> 3     1         3 -483.9464
#> 4     1         4 -482.8378
#> 5     1         5 -484.0765
#> 6     1         6 -484.2357

np <- nuts_params(fit)
str(np)
#> 'data.frame':    24000 obs. of  4 variables:
#>  $ Chain    : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ Iteration: int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ Parameter: Factor w/ 6 levels "accept_stat__",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Value    : num  0.999 0.955 0.964 0.997 0.93 ...
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value)
#> [1] 0

head(rhat(fit))
#>     b_Intercept           b_age     b_sexfemale b_age:sexfemale           sigma 
#>        1.002666        1.002347        1.002364        1.002222        1.001104 
#>       Intercept 
#>        1.001395 
head(neff_ratio(fit))
#>     b_Intercept           b_age     b_sexfemale b_age:sexfemale           sigma 
#>       0.3625533       0.3462164       0.3499752       0.3245068       0.5271566 
#>       Intercept 
#>       0.5210288 
# }
```
