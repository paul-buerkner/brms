# Add model fit criteria to model objects

Add model fit criteria to model objects

## Usage

``` r
add_criterion(x, ...)

# S3 method for class 'brmsfit'
add_criterion(
  x,
  criterion,
  model_name = NULL,
  overwrite = FALSE,
  file = NULL,
  force_save = FALSE,
  ...
)
```

## Arguments

- x:

  An R object typically of class `brmsfit`.

- ...:

  Further arguments passed to the underlying functions computing the
  model fit criteria. If you are recomputing an already stored criterion
  with other `...` arguments, make sure to set `overwrite = TRUE`.

- criterion:

  Names of model fit criteria to compute. Currently supported are
  `"loo"`, `"waic"`, `"kfold"`, `"loo_subsample"`, `"bayes_R2"`
  (Bayesian R-squared), `"loo_R2"` (LOO-adjusted R-squared), and
  `"marglik"` (log marginal likelihood).

- model_name:

  Optional name of the model. If `NULL` (the default) the name is taken
  from the call to `x`.

- overwrite:

  Logical; Indicates if already stored fit indices should be
  overwritten. Defaults to `FALSE`. Setting it to `TRUE` is useful for
  example when changing additional arguments of an already stored
  criterion.

- file:

  Either `NULL` or a character string. In the latter case, the fitted
  model object including the newly added criterion values is saved via
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html) in a file named after
  the string supplied in `file`. The `.rds` extension is added
  automatically. If `x` was already stored in a file before, the file
  name will be reused automatically (with a message) unless overwritten
  by `file`. In any case, `file` only applies if new criteria were
  actually added via `add_criterion` or if `force_save` was set to
  `TRUE`.

- force_save:

  Logical; only relevant if `file` is specified and ignored otherwise.
  If `TRUE`, the fitted model object will be saved regardless of whether
  new criteria were added via `add_criterion`.

## Value

An object of the same class as `x`, but with model fit criteria added
for later usage.

## Details

Functions `add_loo` and `add_waic` are aliases of `add_criterion` with
fixed values for the `criterion` argument.

## Examples

``` r
# \dontrun{
fit <- brm(count ~ Trt, data = epilepsy)
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
#> Chain 1:  Elapsed Time: 0.031 seconds (Warm-up)
#> Chain 1:                0.018 seconds (Sampling)
#> Chain 1:                0.049 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.028 seconds (Warm-up)
#> Chain 2:                0.016 seconds (Sampling)
#> Chain 2:                0.044 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.025 seconds (Warm-up)
#> Chain 3:                0.015 seconds (Sampling)
#> Chain 3:                0.04 seconds (Total)
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
#> Chain 4:  Elapsed Time: 0.028 seconds (Warm-up)
#> Chain 4:                0.016 seconds (Sampling)
#> Chain 4:                0.044 seconds (Total)
#> Chain 4: 
# add both LOO and WAIC at once
fit <- add_criterion(fit, c("loo", "waic"))
#> Warning: Found 1 observations with a pareto_k > 0.7 in model 'fit'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> Warning: 
#> 5 (2.1%) p_waic estimates greater than 0.4. We recommend trying loo instead.
print(fit$criteria$loo)
#> 
#> Computed from 4000 by 236 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -936.1 42.0
#> p_loo        13.5  7.3
#> looic      1872.1 84.1
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.0]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     235   99.6%   338     
#>    (0.7, 1]   (bad)        1    0.4%   <NA>    
#>    (1, Inf)   (very bad)   0    0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.
print(fit$criteria$waic)
#> 
#> Computed from 4000 by 236 log-likelihood matrix.
#> 
#>           Estimate   SE
#> elpd_waic   -936.7 42.5
#> p_waic        14.1  7.9
#> waic        1873.5 85.1
#> 
#> 5 (2.1%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
# }
```
