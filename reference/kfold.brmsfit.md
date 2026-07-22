# K-Fold Cross-Validation

Perform exact K-fold cross-validation by refitting the model \\K\\ times
each leaving out one-\\K\\th of the original data. Folds can be run in
parallel using the future package.

## Usage

``` r
# S3 method for class 'brmsfit'
kfold(
  x,
  ...,
  K = 10,
  Ksub = NULL,
  folds = NULL,
  group = NULL,
  joint = FALSE,
  compare = TRUE,
  resp = NULL,
  model_names = NULL,
  save_fits = FALSE,
  recompile = NULL,
  future_args = list()
)
```

## Arguments

- x:

  A `brmsfit` object.

- ...:

  Further arguments passed to
  [`brm`](https://paulbuerkner.com/brms/reference/brm.md) and
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md).

- K:

  The number of subsets of equal (if possible) size into which the data
  will be partitioned for performing \\K\\-fold cross-validation. The
  model is refit `K` times, each time leaving out one of the `K`
  subsets. If `K` is equal to the total number of observations in the
  data then \\K\\-fold cross-validation is equivalent to exact
  leave-one-out cross-validation.

- Ksub:

  Optional number of subsets (of those subsets defined by `K`) to be
  evaluated. If `NULL` (the default), \\K\\-fold cross-validation will
  be performed on all subsets. If `Ksub` is a single integer, `Ksub`
  subsets (out of all `K`) subsets will be randomly chosen. If `Ksub`
  consists of multiple integers or a one-dimensional array (created via
  `as.array`) potentially of length one, the corresponding subsets will
  be used. This argument is primarily useful, if evaluation of all
  subsets is infeasible for some reason.

- folds:

  Determines how the subsets are being constructed. Possible values are
  `NULL` (the default), `"stratified"`, `"grouped"`, or `"loo"`. May
  also be a vector of length equal to the number of observations in the
  data. Alters the way `group` is handled. More information is provided
  in the 'Details' section.

- group:

  Optional name of a grouping variable or factor in the model. What
  exactly is done with this variable depends on argument `folds`. More
  information is provided in the 'Details' section.

- joint:

  Indicates which observations' log likelihoods shall be considered
  jointly in the ELPD computation. If `"obs"` or `FALSE` (the default),
  each observation is considered separately. This enables comparability
  of `kfold` with `loo`. If `"fold"` or `TRUE`, the joint log
  likelihoods per fold are used. If `"group"`, the joint log likelihoods
  per group within folds are used (only available if argument `group` is
  specified).

- compare:

  A flag indicating if the information criteria of the models should be
  compared to each other via
  [`loo_compare`](https://paulbuerkner.com/brms/reference/loo_compare.brmsfit.md).

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- model_names:

  If `NULL` (the default) will use model names derived from deparsing
  the call. Otherwise will use the passed values as model names.

- save_fits:

  If `TRUE`, a component `fits` is added to the returned object to store
  the cross-validated `brmsfit` objects and the indices of the omitted
  observations for each fold. Defaults to `FALSE`.

- recompile:

  Logical, indicating whether the Stan model should be recompiled. This
  may be necessary if you are running `reloo` on another machine than
  the one used to fit the model.

- future_args:

  A list of further arguments passed to
  [`future`](https://future.futureverse.org/reference/future.html) for
  additional control over parallel execution if activated.

## Value

`kfold` returns an object that has a similar structure as the objects
returned by the `loo` and `waic` methods and can be used with the same
post-processing functions.

## Details

The `kfold` function performs exact \\K\\-fold cross-validation. First
the data are partitioned into \\K\\ folds (i.e. subsets) of equal (or as
close to equal as possible) size by default. Then the model is refit
\\K\\ times, each time leaving out one of the `K` subsets. If \\K\\ is
equal to the total number of observations in the data then \\K\\-fold
cross-validation is equivalent to exact leave-one-out cross-validation
(to which `loo` is an efficient approximation). The `compare_ic`
function is also compatible with the objects returned by `kfold`.

The subsets can be constructed in multiple different ways:

- If both `folds` and `group` are `NULL`, the subsets are randomly
  chosen so that they have equal (or as close to equal as possible)
  size.

- If `folds` is `NULL` but `group` is specified, the data is split up
  into subsets, each time omitting all observations of one of the factor
  levels, while ignoring argument `K`.

- If `folds = "stratified"` the subsets are stratified after `group`
  using
  [`loo::kfold_split_stratified`](https://mc-stan.org/loo/reference/kfold-helpers.html).

- If `folds = "grouped"` the subsets are split by `group` using
  [`loo::kfold_split_grouped`](https://mc-stan.org/loo/reference/kfold-helpers.html).

- If `folds = "loo"` exact leave-one-out cross-validation will be
  performed and `K` will be ignored. Further, if `group` is specified,
  all observations corresponding to the factor level of the currently
  predicted single value are omitted. Thus, in this case, the predicted
  values are only a subset of the omitted ones.

- If `folds` is a numeric vector, it must contain one element per
  observation in the data. Each element of the vector is an integer in
  `1:K` indicating to which of the `K` folds the corresponding
  observation belongs. There are some convenience functions available in
  the loo package that create integer vectors to use for this purpose
  (see the Examples section below and also the
  [kfold-helpers](https://mc-stan.org/loo/reference/kfold-helpers.html)
  page).

When running `kfold` on a `brmsfit` created with the cmdstanr backend in
a different R session, several recompilations will be triggered because
by default, cmdstanr writes the model executable to a temporary
directory. To avoid that, set option `"cmdstanr_write_stan_file_dir"` to
a nontemporary path of your choice before creating the original
`brmsfit` (see section 'Examples' below).

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
[`reloo`](https://paulbuerkner.com/brms/reference/reloo.brmsfit.md)

## Examples

``` r
# \dontrun{
fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
           data = epilepsy, family = poisson())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.48 seconds.
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
#> Chain 1:  Elapsed Time: 4.27 seconds (Warm-up)
#> Chain 1:                2.574 seconds (Sampling)
#> Chain 1:                6.844 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.41 seconds.
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
#> Chain 2:  Elapsed Time: 4.693 seconds (Warm-up)
#> Chain 2:                2.522 seconds (Sampling)
#> Chain 2:                7.215 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4.1e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.41 seconds.
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
#> Chain 3:  Elapsed Time: 4.672 seconds (Warm-up)
#> Chain 3:                4.459 seconds (Sampling)
#> Chain 3:                9.131 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.4 seconds.
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
#> Chain 4:  Elapsed Time: 4.52 seconds (Warm-up)
#> Chain 4:                2.463 seconds (Sampling)
#> Chain 4:                6.983 seconds (Total)
#> Chain 4: 
# throws warning about some pareto k estimates being too high
(loo1 <- loo(fit1))
#> Warning: Found 56 observations with a pareto_k > 0.7 in model 'fit1'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> 
#> Computed from 4000 by 236 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -593.9 13.7
#> p_loo       106.6  6.9
#> looic      1187.9 27.4
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.7]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     180   76.3%   195     
#>    (0.7, 1]   (bad)       52   22.0%   <NA>    
#>    (1, Inf)   (very bad)   4    1.7%   <NA>    
#> See help('pareto-k-diagnostic') for details.
# perform 10-fold cross validation
(kfold1 <- kfold(fit1, chains = 1))
#> Fitting model 1 out of 10
#> Start sampling
#> Fitting model 2 out of 10
#> Start sampling
#> Fitting model 3 out of 10
#> Start sampling
#> Fitting model 4 out of 10
#> Start sampling
#> Fitting model 5 out of 10
#> Start sampling
#> Fitting model 6 out of 10
#> Start sampling
#> Fitting model 7 out of 10
#> Start sampling
#> Fitting model 8 out of 10
#> Start sampling
#> Fitting model 9 out of 10
#> Start sampling
#> Fitting model 10 out of 10
#> Start sampling
#> Warning: Found 3 observations with a pareto_k > 0.666666666666667 in model 'fit1'.
#> 
#> Based on 10-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -618.4 16.9
#> p_kfold       131.0 10.4
#> kfoldic      1236.7 33.9
#> ------
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.67]   (good)     233   98.7%   556     
#>    (0.67, 1]   (bad)        0    0.0%   <NA>    
#>     (1, Inf)   (very bad)   3    1.3%   <NA>    
#> See help('pareto-k-diagnostic') for details.

# use joint likelihoods per fold for ELPD evaluation
kfold(fit1, chains = 1, joint = "fold")
#> Fitting model 1 out of 10
#> Start sampling
#> Fitting model 2 out of 10
#> Start sampling
#> Fitting model 3 out of 10
#> Start sampling
#> Fitting model 4 out of 10
#> Start sampling
#> Fitting model 5 out of 10
#> Start sampling
#> Fitting model 6 out of 10
#> Start sampling
#> Fitting model 7 out of 10
#> Start sampling
#> Fitting model 8 out of 10
#> Start sampling
#> Fitting model 9 out of 10
#> Start sampling
#> Fitting model 10 out of 10
#> Start sampling
#> Warning: Found 10 observations with a pareto_k > 0.666666666666667 in model 'fit1'.
#> 
#> Based on 10-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -656.0 22.4
#> p_kfold       169.3 18.1
#> kfoldic      1312.0 44.8
#> ------
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.67]   (good)      0      0.0%  <NA>    
#>    (0.67, 1]   (bad)       0      0.0%  <NA>    
#>     (1, Inf)   (very bad) 10    100.0%  <NA>    
#> See help('pareto-k-diagnostic') for details.

# use the future package for parallelization of models
# that is to fit models belonging to different folds in parallel
library(future)
plan(multisession, workers = 4)
kfold(fit1, chains = 1)
#> Fitting model 1 out of 10
#> Start sampling
#> Fitting model 2 out of 10
#> Start sampling
#> Fitting model 3 out of 10
#> Start sampling
#> Fitting model 4 out of 10
#> Start sampling
#> Fitting model 5 out of 10
#> Start sampling
#> Fitting model 6 out of 10
#> Start sampling
#> Fitting model 7 out of 10
#> Start sampling
#> Fitting model 8 out of 10
#> Start sampling
#> Fitting model 9 out of 10
#> Start sampling
#> Fitting model 10 out of 10
#> Start sampling
#> Warning: Found 2 observations with a pareto_k > 0.666666666666667 in model 'fit1'.
#> 
#> Based on 10-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -616.1 16.8
#> p_kfold       128.8 10.2
#> kfoldic      1232.2 33.6
#> ------
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.67]   (good)     234   99.2%   599     
#>    (0.67, 1]   (bad)        0    0.0%   <NA>    
#>     (1, Inf)   (very bad)   2    0.8%   <NA>    
#> See help('pareto-k-diagnostic') for details.
plan(sequential)

## to avoid recompilations when running kfold() on a 'cmdstanr'-backend fit
## in a fresh R session, set option 'cmdstanr_write_stan_file_dir' before
## creating the initial 'brmsfit'
## CAUTION: the following code creates some files in the current working
## directory: two 'model_<hash>.stan' files, one 'model_<hash>(.exe)'
## executable, and one 'fit_cmdstanr_<some_number>.rds' file
set.seed(7)
fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
options(cmdstanr_write_stan_file_dir = getwd())
fit_cmdstanr <- brm(rate ~ conc + state, data = Puromycin,
                    backend = "cmdstanr", file = fname)
#> Start sampling
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 1 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 1 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 1 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 1 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 1 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 1 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 1 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 1 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 1 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 1 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 2 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 2 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 2 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 2 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 2 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 2 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 2 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 2 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 2 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 2 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 2 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 2 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 2 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 2 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 2 finished in 0.0 seconds.
#> Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 3 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 3 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 3 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 3 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 3 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 3 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 3 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 3 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 3 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 3 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 3 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 3 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 3 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 3 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 3 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 3 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 4 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 4 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 4 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 4 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 4 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 4 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 4 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 4 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 4 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 4 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 4 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 4 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 4 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 4 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 4 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 4 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.9 seconds.
#> 

# now restart the R session and run the following (after attaching 'brms')
set.seed(7)
fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
fit_cmdstanr <- brm(rate ~ conc + state,
                    data = Puromycin,
                    backend = "cmdstanr",
                    file = fname)
kfold_cmdstanr <- kfold(fit_cmdstanr, K = 2)
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 finished in 0.0 seconds.
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.8 seconds.
#> 
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 finished in 0.1 seconds.
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.8 seconds.
#> 
#> Fitting model 1 out of 2
#> Start sampling
#> Fitting model 2 out of 2
#> Start sampling
# }
```
