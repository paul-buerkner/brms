# Pareto smoothed importance sampling (PSIS)

Implementation of Pareto smoothed importance sampling (PSIS), a method
for stabilizing importance ratios. The version of PSIS implemented here
corresponds to the algorithm presented in Vehtari, Simpson, Gelman, Yao,
and Gabry (2024). For PSIS diagnostics see the
[pareto-k-diagnostic](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
page.

## Usage

``` r
# S3 method for class 'brmsfit'
psis(log_ratios, newdata = NULL, resp = NULL, model_name = NULL, ...)
```

## Arguments

- log_ratios:

  A fitted model object of class `brmsfit`. Argument is named
  "log_ratios" to match the argument name of the
  [`loo::psis`](https://mc-stan.org/loo/reference/psis.html) generic
  function.

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

- model_name:

  Currently ignored.

- ...:

  Further arguments passed to
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md)
  and [`loo::psis`](https://mc-stan.org/loo/reference/psis.html).

## Value

The `psis()` methods return an object of class `"psis"`, which is a
named list with the following components:

- `log_weights`:

  Vector or matrix of smoothed (and truncated) but *unnormalized* log
  weights. To get normalized weights use the
  [`weights()`](https://mc-stan.org/loo/reference/weights.importance_sampling.html)
  method provided for objects of class `"psis"`.

- `diagnostics`:

  A named list containing two vectors:

  - `pareto_k`: Estimates of the shape parameter \\k\\ of the
    generalized Pareto distribution. See the
    [pareto-k-diagnostic](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
    page for details.

  - `n_eff`: PSIS effective sample size estimates.

Objects of class `"psis"` also have the following
[attributes](https://rdrr.io/r/base/attributes.html):

- `norm_const_log`:

  Vector of precomputed values of `colLogSumExps(log_weights)` that are
  used internally by the `weights` method to normalize the log weights.

- `tail_len`:

  Vector of tail lengths used for fitting the generalized Pareto
  distribution.

- `r_eff`:

  If specified, the user's `r_eff` argument.

- `dims`:

  Integer vector of length 2 containing `S` (posterior sample size) and
  `N` (number of observations).

- `method`:

  Method used for importance sampling, here `psis`.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4
([journal
version](https://link.springer.com/article/10.1007/s11222-016-9696-4),
[preprint arXiv:1507.04544](https://arxiv.org/abs/1507.04544)).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72):1-58. [PDF](https://jmlr.org/papers/v25/19-556.html)

## Examples

``` r
# \dontrun{
fit <- brm(rating ~ treat + period + carry, data = inhaler)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.14 seconds.
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
#> Chain 1:  Elapsed Time: 0.03 seconds (Warm-up)
#> Chain 1:                0.028 seconds (Sampling)
#> Chain 1:                0.058 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 6e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
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
#> Chain 2:  Elapsed Time: 0.03 seconds (Warm-up)
#> Chain 2:                0.027 seconds (Sampling)
#> Chain 2:                0.057 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 6e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
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
#> Chain 3:  Elapsed Time: 0.03 seconds (Warm-up)
#> Chain 3:                0.029 seconds (Sampling)
#> Chain 3:                0.059 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 6e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
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
#> Chain 4:  Elapsed Time: 0.029 seconds (Warm-up)
#> Chain 4:                0.028 seconds (Sampling)
#> Chain 4:                0.057 seconds (Total)
#> Chain 4: 
psis(fit)
#> Computed from 4000 by 572 log-weights matrix.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [1.0, 1.2]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
# }
```
