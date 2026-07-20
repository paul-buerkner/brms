# Efficient approximate leave-one-out cross-validation (LOO)

Perform approximate leave-one-out cross-validation based on the
posterior likelihood using the loo package. For more details see
[`loo`](https://mc-stan.org/loo/reference/loo.html).

## Usage

``` r
# S3 method for class 'brmsfit'
loo(
  x,
  ...,
  compare = TRUE,
  resp = NULL,
  pointwise = FALSE,
  moment_match = FALSE,
  reloo = FALSE,
  k_threshold = 0.7,
  save_psis = FALSE,
  moment_match_args = list(),
  reloo_args = list(),
  model_names = NULL
)
```

## Arguments

- x:

  A `brmsfit` object.

- ...:

  More `brmsfit` objects or further arguments passed to the underlying
  post-processing functions. In particular, see
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  for further supported arguments.

- compare:

  A flag indicating if the information criteria of the models should be
  compared to each other via
  [`loo_compare`](https://paulbuerkner.com/brms/reference/loo_compare.brmsfit.md).

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- pointwise:

  A flag indicating whether to compute the full log-likelihood matrix at
  once or separately for each observation. The latter approach is
  usually considerably slower but requires much less working memory.
  Accordingly, if one runs into memory issues, `pointwise = TRUE` is the
  way to go.

- moment_match:

  Logical; Indicate whether
  [`loo_moment_match`](https://paulbuerkner.com/brms/reference/loo_moment_match.brmsfit.md)
  should be applied on problematic observations. Defaults to `FALSE`.
  For most models, moment matching will only work if you have set
  `save_pars = save_pars(all = TRUE)` when fitting the model with
  [`brm`](https://paulbuerkner.com/brms/reference/brm.md). See
  [`loo_moment_match.brmsfit`](https://paulbuerkner.com/brms/reference/loo_moment_match.brmsfit.md)
  for more details.

- reloo:

  Logical; Indicate whether
  [`reloo`](https://paulbuerkner.com/brms/reference/reloo.brmsfit.md)
  should be applied on problematic observations. Defaults to `FALSE`.

- k_threshold:

  The Pareto \\k\\ threshold for which observations
  [`loo_moment_match`](https://paulbuerkner.com/brms/reference/loo_moment_match.brmsfit.md)
  or [`reloo`](https://paulbuerkner.com/brms/reference/reloo.brmsfit.md)
  is applied if argument `moment_match` or `reloo` is `TRUE`. Defaults
  to `0.7`. See
  [`pareto_k_ids`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
  for more details.

- save_psis:

  Should the `"psis"` object created internally be saved in the returned
  object? For more details see
  [`loo`](https://mc-stan.org/loo/reference/loo.html).

- moment_match_args:

  Optional named `list` of additional arguments passed to
  [`loo_moment_match`](https://paulbuerkner.com/brms/reference/loo_moment_match.brmsfit.md).

- reloo_args:

  Optional named `list` of additional arguments passed to
  [`reloo`](https://paulbuerkner.com/brms/reference/reloo.brmsfit.md).
  This can be useful, among others, to control how many chains,
  iterations, etc. to use for the fitted sub-models.

- model_names:

  If `NULL` (the default) will use model names derived from deparsing
  the call. Otherwise will use the passed values as model names.

## Value

If just one object is provided, an object of class `loo`. If multiple
objects are provided, an object of class `loolist`.

## Details

See
[`loo_compare`](https://paulbuerkner.com/brms/reference/loo_compare.brmsfit.md)
for details on model comparisons. For `brmsfit` objects, `LOO` is an
alias of `loo`. Use method
[`add_criterion`](https://paulbuerkner.com/brms/reference/add_criterion.md)
to store information criteria in the fitted model object for later
usage.

## References

Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. In Statistics
and Computing, doi:10.1007/s11222-016-9696-4. arXiv preprint
arXiv:1507.04544.

Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models. Statistics and Computing, 24,
997-1016.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
The Journal of Machine Learning Research, 11, 3571-3594.

## Examples

``` r
# \dontrun{
# model with population-level effects only
fit1 <- brm(rating ~ treat + period + carry,
            data = inhaler)
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
#> Chain 1:  Elapsed Time: 0.035 seconds (Warm-up)
#> Chain 1:                0.032 seconds (Sampling)
#> Chain 1:                0.067 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 5e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
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
#> Chain 2:  Elapsed Time: 0.035 seconds (Warm-up)
#> Chain 2:                0.032 seconds (Sampling)
#> Chain 2:                0.067 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 5e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
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
#> Chain 3:  Elapsed Time: 0.035 seconds (Warm-up)
#> Chain 3:                0.032 seconds (Sampling)
#> Chain 3:                0.067 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 5e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
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
#> Chain 4:  Elapsed Time: 0.035 seconds (Warm-up)
#> Chain 4:                0.032 seconds (Sampling)
#> Chain 4:                0.067 seconds (Total)
#> Chain 4: 
(loo1 <- loo(fit1))
#> 
#> Computed from 4000 by 572 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -529.5 25.9
#> p_loo         6.3  1.0
#> looic      1059.0 51.9
#> ------
#> MCSE of elpd_loo is 0.0.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [1.0, 1.3]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.

# model with an additional varying intercept for subjects
fit2 <- brm(rating ~ treat + period + carry + (1|subject),
            data = inhaler)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.54 seconds.
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
#> Chain 1:  Elapsed Time: 1.645 seconds (Warm-up)
#> Chain 1:                0.778 seconds (Sampling)
#> Chain 1:                2.423 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4.5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 2:  Elapsed Time: 1.704 seconds (Warm-up)
#> Chain 2:                0.782 seconds (Sampling)
#> Chain 2:                2.486 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4.5e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 3:  Elapsed Time: 1.568 seconds (Warm-up)
#> Chain 3:                0.777 seconds (Sampling)
#> Chain 3:                2.345 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4.5e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 4:  Elapsed Time: 1.582 seconds (Warm-up)
#> Chain 4:                0.778 seconds (Sampling)
#> Chain 4:                2.36 seconds (Total)
#> Chain 4: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
(loo2 <- loo(fit2))
#> Warning: Found 4 observations with a pareto_k > 0.7 in model 'fit2'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> 
#> Computed from 4000 by 572 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -522.0 26.3
#> p_loo        87.3  7.7
#> looic      1044.0 52.6
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 1.9]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     568   99.3%   106     
#>    (0.7, 1]   (bad)        4    0.7%   <NA>    
#>    (1, Inf)   (very bad)   0    0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.

# compare both models
loo_compare(loo1, loo2)
#>  model elpd_diff se_diff p_worse diag_diff      diag_elpd
#>   fit2       0.0     0.0      NA           4 k_psis > 0.7
#>   fit1      -7.5     4.7    0.95                         
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
