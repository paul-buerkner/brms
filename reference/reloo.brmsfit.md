# Compute exact cross-validation for problematic observations

Compute exact cross-validation for problematic observations for which
approximate leave-one-out cross-validation may return incorrect results.
Models for problematic observations can be run in parallel using the
future package.

## Usage

``` r
# S3 method for class 'brmsfit'
reloo(
  x,
  loo = NULL,
  k_threshold = 0.7,
  newdata = NULL,
  resp = NULL,
  check = TRUE,
  recompile = NULL,
  future_args = list(),
  ...
)

# S3 method for class 'loo'
reloo(x, fit, ...)

reloo(x, ...)
```

## Arguments

- x:

  An R object of class `brmsfit` or `loo` depending on the method.

- loo:

  An R object of class `loo`. If `NULL`, brms will try to extract a
  precomputed `loo` object from the fitted model, added there via
  [`add_criterion`](https://paulbuerkner.com/brms/reference/add_criterion.md).

- k_threshold:

  The threshold at which Pareto \\k\\ estimates are treated as
  problematic. Defaults to `0.7`. See
  [`pareto_k_ids`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
  for more details.

- newdata:

  An optional data.frame for which to evaluate predictions. If `NULL`
  (default), the original data of the model is used. `NA` values within
  factors (excluding grouping variables) are interpreted as if all dummy
  variables of this factor are zero. This allows, for instance, to make
  predictions of the grand mean when using sum coding. `NA` values
  within grouping variables are treated as a new level.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- check:

  Logical; If `TRUE` (the default), some checks check are performed if
  the `loo` object was generated from the `brmsfit` object passed to
  argument `fit`.

- recompile:

  Logical, indicating whether the Stan model should be recompiled. This
  may be necessary if you are running `reloo` on another machine than
  the one used to fit the model.

- future_args:

  A list of further arguments passed to
  [`future`](https://future.futureverse.org/reference/future.html) for
  additional control over parallel execution if activated.

- ...:

  Further arguments passed to
  [`update.brmsfit`](https://paulbuerkner.com/brms/reference/update.brmsfit.md)
  and
  [`log_lik.brmsfit`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md).

- fit:

  An R object of class `brmsfit`.

## Value

An object of the class `loo`.

## Details

Warnings about Pareto \\k\\ estimates indicate observations for which
the approximation to LOO is problematic (this is described in detail in
Vehtari, Gelman, and Gabry (2017) and the
[loo](https://mc-stan.org/loo/reference/loo-package.html) package
documentation). If there are \\J\\ observations with \\k\\ estimates
above `k_threshold`, then `reloo` will refit the original model \\J\\
times, each time leaving out one of the \\J\\ problematic observations.
The pointwise contributions of these observations to the total ELPD are
then computed directly and substituted for the previous estimates from
these \\J\\ observations that are stored in the original `loo` object.

By default, this method uses `sample_new_levels = "gaussian"` to sample
parameter values for new grouping-factor levels (see also
[`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)).
This default will fail for models with non-Gaussian group-level effects.
In this case, we recommend setting `sample_new_levels = "uncertainty"`.

## Parallelization with multiple CPU cores

brms can make use of multiple CPU cores in parallel to speed up
computations in various ways. For efficient use of the available
resources it is recommended to only use parallelism to an extend such
that the available physical CPUs are not oversubscribed. For example,
when you have 8 CPU cores locally available, then you may consider to
run 4 chains with 2 threads per chain for best performance if you happen
to just run a single model. In case you run a simulation study which
requires to run many times a given model, then neither chain nor
within-chain parallelization is advisable as the computational resources
are already exhausted by the simulation study and any further
parallelization beyond the simulation study itself will in fact slow
down the overall runtime. Please be aware that for historical reasons
the nomenclature of the arguments is possibly confusing. The `cores`
argument refers to running different chains in parallel and the
within-chain parallelization will allocate for each chain as many
threads as requested. The requested threads therefore increase the use
of overall CPUs in a multiplicative way.

For more advanced parallelization (including beyond single model fits),
brms also integrates with the future package. Importantly, this enables
seamless integration with the mirai parallelization framework through
the use of the future.mirai adapter. With mirai local and remote
machines can be used in a fully transparent manner to the user. This
includes the possibility to use large number of remote machines running
in the context of a computer cluster, which are managed with queuing
systems. Please refer to the section on distributed computing of
[`mirai::daemons`](https://mirai.r-lib.org/reference/daemons.html).

## See also

[`loo`](https://paulbuerkner.com/brms/reference/loo.brmsfit.md),
[`kfold`](https://paulbuerkner.com/brms/reference/kfold.brmsfit.md)

## Examples

``` r
# \dontrun{
fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
            data = epilepsy, family = poisson())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.57 seconds.
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
#> Chain 1:  Elapsed Time: 2.077 seconds (Warm-up)
#> Chain 1:                1.904 seconds (Sampling)
#> Chain 1:                3.981 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.6e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
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
#> Chain 2:  Elapsed Time: 1.902 seconds (Warm-up)
#> Chain 2:                1.426 seconds (Sampling)
#> Chain 2:                3.328 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.9e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.29 seconds.
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
#> Chain 3:  Elapsed Time: 2.028 seconds (Warm-up)
#> Chain 3:                1.63 seconds (Sampling)
#> Chain 3:                3.658 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.3 seconds.
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
#> Chain 4:  Elapsed Time: 1.965 seconds (Warm-up)
#> Chain 4:                1.428 seconds (Sampling)
#> Chain 4:                3.393 seconds (Total)
#> Chain 4: 

# throws warning about some pareto k estimates being too high
(loo1 <- loo(fit1))
#> Warning: Found 7 observations with a pareto_k > 0.7 in model 'fit1'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> 
#> Computed from 4000 by 236 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -670.4 36.0
#> p_loo        93.3 14.2
#> looic      1340.8 72.1
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.3]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     229   97.0%   164     
#>    (0.7, 1]   (bad)        6    2.5%   <NA>    
#>    (1, Inf)   (very bad)   1    0.4%   <NA>    
#> See help('pareto-k-diagnostic') for details.

# no more warnings after reloo
(reloo1 <- reloo(fit1, loo = loo1, chains = 1))
#> 7 problematic observation(s) found.
#> The model will be refit 7 times.
#> 
#> Fitting model 1 out of 7 (leaving out observation 8)
#> Start sampling
#> 
#> Fitting model 2 out of 7 (leaving out observation 10)
#> Start sampling
#> 
#> Fitting model 3 out of 7 (leaving out observation 16)
#> Start sampling
#> 
#> Fitting model 4 out of 7 (leaving out observation 74)
#> Start sampling
#> 
#> Fitting model 5 out of 7 (leaving out observation 98)
#> Start sampling
#> 
#> Fitting model 6 out of 7 (leaving out observation 115)
#> Start sampling
#> 
#> Fitting model 7 out of 7 (leaving out observation 143)
#> Start sampling
#> 
#> Computed from 4000 by 236 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -672.8 37.3
#> p_loo        95.8 15.6
#> looic      1345.7 74.6
#> ------
#> MCSE of elpd_loo is 0.6.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.3]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
# }
```
