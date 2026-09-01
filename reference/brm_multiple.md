# Run the same brms model on multiple datasets

Run the same brms model on multiple datasets and then combine the
results into one fitted model object. This is useful in particular for
multiple missing value imputation, where the same model is fitted on
multiple imputed data sets. Models can be run in parallel using the
future package.

## Usage

``` r
brm_multiple(
  formula,
  data,
  family = gaussian(),
  prior = NULL,
  data2 = NULL,
  autocor = NULL,
  cov_ranef = NULL,
  sample_prior = c("no", "yes", "only"),
  sparse = NULL,
  knots = NULL,
  stanvars = NULL,
  stan_funs = NULL,
  silent = getOption("brms.silent", 1),
  recompile = FALSE,
  combine = TRUE,
  fit = NA,
  algorithm = getOption("brms.algorithm", "sampling"),
  seed = NA,
  file = NULL,
  file_compress = TRUE,
  file_refit = getOption("brms.file_refit", "never"),
  ...
)
```

## Arguments

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html),
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
  or
  [`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)
  (or one that can be coerced to that classes): A symbolic description
  of the model to be fitted. The details of model specification are
  explained in
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

- data:

  A *list* of data.frames each of which will be used to fit a separate
  model. Alternatively, a `mids` object from the mice package.

- family:

  A description of the response distribution and link function to be
  used in the model. This can be a family function, a call to a family
  function or a character string naming the family. Every family
  function has a `link` argument allowing to specify the link function
  to be applied on the response variable. If not specified, default
  links are used. For details of supported families see
  [`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).
  By default, a linear `gaussian` model is applied. In multivariate
  models, `family` might also be a list of families.

- prior:

  One or more `brmsprior` objects created by
  [`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md) or
  related functions and combined using the `c` method or the `+`
  operator. See also
  [`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md)
  for more help.

- data2:

  A *list* of named lists each of which will be used to fit a separate
  model. Each of the named lists contains objects representing data
  which cannot be passed via argument `data` (see
  [`brm`](https://paulbuerkner.com/brms/reference/brm.md) for examples).
  The length of the outer list should match the length of the list
  passed to the `data` argument.

- autocor:

  (Deprecated) An optional
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md)
  object describing the correlation structure within the response
  variable (i.e., the 'autocorrelation'). See the documentation of
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md) for
  a description of the available correlation structures. Defaults to
  `NULL`, corresponding to no correlations. In multivariate models,
  `autocor` might also be a list of autocorrelation structures. It is
  now recommend to specify autocorrelation terms directly within
  `formula`. See
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  for more details.

- cov_ranef:

  (Deprecated) A list of matrices that are proportional to the (within)
  covariance structure of the group-level effects. The names of the
  matrices should correspond to columns in `data` that are used as
  grouping factors. All levels of the grouping factor should appear as
  rownames of the corresponding matrix. This argument can be used, among
  others to model pedigrees and phylogenetic effects. It is now
  recommended to specify those matrices in the formula interface using
  the [`gr`](https://paulbuerkner.com/brms/reference/gr.md) and related
  functions. See
  [`vignette("brms_phylogenetics")`](https://paulbuerkner.com/brms/articles/brms_phylogenetics.md)
  for more details.

- sample_prior:

  Indicate if draws from priors should be drawn additionally to the
  posterior draws. Options are `"no"` (the default), `"yes"`, and
  `"only"`. Among others, these draws can be used to calculate Bayes
  factors for point hypotheses via
  [`hypothesis`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.md).
  Please note that improper priors are not sampled, including the
  default improper priors used by `brm`. See
  [`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md) on
  how to set (proper) priors. Please also note that prior draws for the
  overall intercept are not obtained by default for technical reasons.
  See
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  how to obtain prior draws for the intercept. If `sample_prior` is set
  to `"only"`, draws are drawn solely from the priors ignoring the
  likelihood, which allows among others to generate draws from the prior
  predictive distribution. In this case, all parameters must have proper
  priors.

- sparse:

  (Deprecated) Logical; indicates whether the population-level design
  matrices should be treated as sparse (defaults to `FALSE`). For design
  matrices with many zeros, this can considerably reduce required
  memory. Sampling speed is currently not improved or even slightly
  decreased. It is now recommended to use the `sparse` argument of
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  and related functions.

- knots:

  Optional list containing user specified knot values to be used for
  basis construction of smoothing terms. See
  [`gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html) for more details.

- stanvars:

  An optional `stanvars` object generated by function
  [`stanvar`](https://paulbuerkner.com/brms/reference/stanvar.md) to
  define additional variables for use in Stan's program blocks.

- stan_funs:

  (Deprecated) An optional character string containing self-defined Stan
  functions, which will be included in the functions block of the
  generated Stan code. It is now recommended to use the `stanvars`
  argument for this purpose instead.

- silent:

  Verbosity level between `0` and `2`. If `1` (the default), most of the
  informational messages of compiler and sampler are suppressed. If `2`,
  even more messages are suppressed. The actual sampling progress is
  still printed. Set `refresh = 0` to turn this off as well. If using
  `backend = "rstan"` you can also set `open_progress = FALSE` to
  prevent opening additional progress bars. Can be set globally for the
  current R session via the `"brms.silent"` option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- recompile:

  Logical, indicating whether the Stan model should be recompiled for
  every imputed data set. Defaults to `FALSE`. If `NULL`, `brm_multiple`
  tries to figure out internally, if recompilation is necessary, for
  example because data-dependent priors have changed. Using the default
  of no recompilation should be fine in most cases.

- combine:

  Logical; Indicates if the fitted models should be combined into a
  single fitted model object via
  [`combine_models`](https://paulbuerkner.com/brms/reference/combine_models.md).
  Defaults to `TRUE`.

- fit:

  An instance of S3 class `brmsfit_multiple` derived from a previous
  fit; defaults to `NA`. If `fit` is of class `brmsfit_multiple`, the
  compiled model associated with the fitted result is re-used and all
  arguments modifying the model code or data are ignored. It is not
  recommended to use this argument directly, but to call the
  [`update`](https://paulbuerkner.com/brms/reference/update.brmsfit_multiple.md)
  method, instead.

- algorithm:

  Character string naming the estimation approach to use. Options are
  `"sampling"` for MCMC (the default), `"meanfield"` for variational
  inference with independent normal distributions, `"fullrank"` for
  variational inference with a multivariate normal distribution,
  `"pathfinder"` for the pathfinder algorithm, `"laplace"` for the
  laplace approximation, or `"fixed_param"` for sampling from fixed
  parameter values. Can be set globally for the current R session via
  the `"brms.algorithm"` option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- seed:

  The seed for random number generation to make results reproducible. If
  `NA` (the default), Stan will set the seed randomly.

- file:

  Either `NULL` or a character string. In the latter case, the fitted
  model object is saved via
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html) in a file named after
  the string supplied in `file`. The `.rds` extension is added
  automatically. If the file already exists, `brm` will load and return
  the saved model object instead of refitting the model. Unless you
  specify the `file_refit` argument as well, the existing files won't be
  overwritten, you have to manually remove the file in order to refit
  and save the model under an existing file name. The file name is
  stored in the `brmsfit` object for later usage.

- file_compress:

  Logical or a character string, specifying one of the compression
  algorithms supported by
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html). If the `file`
  argument is provided, this compression will be used when saving the
  fitted model object.

- file_refit:

  Modifies when the fit stored via the `file` argument is re-used. Can
  be set globally for the current R session via the `"brms.file_refit"`
  option (see [`options`](https://rdrr.io/r/base/options.html)). For
  `"never"` (default) the fit is always loaded if it exists and fitting
  is skipped. For `"always"` the model is always refitted. If set to
  `"on_change"`, brms will refit the model if model, data or algorithm
  as passed to Stan differ from what is stored in the file. This also
  covers changes in priors, `sample_prior`, `stanvars`, covariance
  structure, etc. If you believe there was a false positive, you can use
  [`brmsfit_needs_refit`](https://paulbuerkner.com/brms/reference/brmsfit_needs_refit.md)
  to see why refit is deemed necessary. Refit will not be triggered for
  changes in additional parameters of the fit (e.g., initial values,
  number of iterations, control arguments, ...). A known limitation is
  that a refit will be triggered if within-chain parallelization is
  switched on/off.

- ...:

  Further arguments passed to
  [`brm`](https://paulbuerkner.com/brms/reference/brm.md).

## Value

If `combine = TRUE` a `brmsfit_multiple` object, which inherits from
class `brmsfit` and behaves essentially the same. If `combine = FALSE` a
list of `brmsfit` objects.

## Details

The inference for the combined model posterior may issue false positive
convergence warnings, as the MCMC chains corresponding to posteriors
with different datasets may not necessarily overlap, even if the
inference for each of the original posterior did converge. To find out
whether the inference for each of the original posterior converged,
subset the draws belonging to the individual posteriors (model fits) and
then run convergence diagnostics. See Examples below for details.

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

## Examples

``` r
# \dontrun{
library(mice)
#> 
#> Attaching package: ‘mice’
#> The following object is masked from ‘package:stats’:
#> 
#>     filter
#> The following objects are masked from ‘package:base’:
#> 
#>     cbind, rbind
m <- 5
imp <- mice(nhanes2, m = m)
#> 
#>  iter imp variable
#>   1   1  bmi  hyp  chl
#>   1   2  bmi  hyp  chl
#>   1   3  bmi  hyp  chl
#>   1   4  bmi  hyp  chl
#>   1   5  bmi  hyp  chl
#>   2   1  bmi  hyp  chl
#>   2   2  bmi  hyp  chl
#>   2   3  bmi  hyp  chl
#>   2   4  bmi  hyp  chl
#>   2   5  bmi  hyp  chl
#>   3   1  bmi  hyp  chl
#>   3   2  bmi  hyp  chl
#>   3   3  bmi  hyp  chl
#>   3   4  bmi  hyp  chl
#>   3   5  bmi  hyp  chl
#>   4   1  bmi  hyp  chl
#>   4   2  bmi  hyp  chl
#>   4   3  bmi  hyp  chl
#>   4   4  bmi  hyp  chl
#>   4   5  bmi  hyp  chl
#>   5   1  bmi  hyp  chl
#>   5   2  bmi  hyp  chl
#>   5   3  bmi  hyp  chl
#>   5   4  bmi  hyp  chl
#>   5   5  bmi  hyp  chl

# fit the model using mice and lm
fit_imp1 <- with(lm(bmi ~ age + hyp + chl), data = imp)
summary(pool(fit_imp1))
#>          term   estimate std.error statistic       df     p.value
#> 1 (Intercept) 19.0050191 3.9623275  4.796428 7.392980 0.001696942
#> 2    age40-59 -5.7819661 2.5790688 -2.241881 3.758055 0.092760029
#> 3    age60-99 -8.0118093 3.3519357 -2.390204 3.173096 0.092008489
#> 4      hypyes  2.2912320 2.2194912  1.032323 6.733643 0.337565925
#> 5         chl  0.0563835 0.0237295  2.376093 5.516864 0.058724384

# fit the model using brms
fit_imp2 <- brm_multiple(bmi ~ age + hyp + chl, data = imp, chains = 1)
#> Compiling the C++ model
#> Error: ‘nbrOfWorkers >= 1L’ is not TRUE
summary(fit_imp2)
#> Error: object 'fit_imp2' not found
plot(fit_imp2, variable = "^b_", regex = TRUE)
#> Error: object 'fit_imp2' not found

# investigate convergence of inference for the original posteriors
library(posterior)
#> This is posterior version 1.7.0
#> 
#> Attaching package: ‘posterior’
#> The following objects are masked from ‘package:stats’:
#> 
#>     mad, sd, var
#> The following objects are masked from ‘package:base’:
#> 
#>     %in%, match
draws <- as_draws_array(fit_imp2)
#> Error: object 'fit_imp2' not found
# every dataset has just one chain here
draws_per_dat <- lapply(1:m, \(i) subset_draws(draws, chain = i))
#> Error in FUN(X[[i]], ...): object 'draws' not found
lapply(draws_per_dat, summarise_draws, default_convergence_measures())
#> Error: object 'draws_per_dat' not found

# use the future package for parallelization
library(future)
plan(multisession, workers = 4)
fit_imp3 <- brm_multiple(bmi ~ age + hyp + chl, data = imp, chains = 1)
#> Compiling the C++ model
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.11 seconds.
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
#> Chain 1:  Elapsed Time: 0.072 seconds (Warm-up)
#> Chain 1:                0.017 seconds (Sampling)
#> Chain 1:                0.089 seconds (Total)
#> Chain 1: 
#> Fitting imputed model 1
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.11 seconds.
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
#> Chain 1:  Elapsed Time: 0.059 seconds (Warm-up)
#> Chain 1:                0.017 seconds (Sampling)
#> Chain 1:                0.076 seconds (Total)
#> Chain 1: 
#> Fitting imputed model 2
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.11 seconds.
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
#> Chain 1:  Elapsed Time: 0.059 seconds (Warm-up)
#> Chain 1:                0.016 seconds (Sampling)
#> Chain 1:                0.075 seconds (Total)
#> Chain 1: 
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
#> Chain 1:  Elapsed Time: 0.065 seconds (Warm-up)
#> Chain 1:                0.018 seconds (Sampling)
#> Chain 1:                0.083 seconds (Total)
#> Chain 1: 
#> Fitting imputed model 3
#> Start sampling
#> Fitting imputed model 4
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.12 seconds.
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
#> Chain 1:  Elapsed Time: 0.09 seconds (Warm-up)
#> Chain 1:                0.028 seconds (Sampling)
#> Chain 1:                0.118 seconds (Total)
#> Chain 1: 
#> Fitting imputed model 5
#> Start sampling
summary(fit_imp3)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: bmi ~ age + hyp + chl 
#>    Data: imp (Number of observations: 25) 
#>   Draws: 5 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 5000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI
#> Intercept    19.06      3.75    11.78    26.20
#> age40M59     -5.81      2.31   -10.22    -1.12
#> age60M99     -8.02      2.92   -13.71    -2.67
#> hypyes        2.32      2.11    -1.98     6.41
#> chl           0.06      0.02     0.01     0.10
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     2.97      0.55     2.06     4.22
#> 
#> Draws were sampled using sampling(NUTS). Overall Rhat and ESS estimates
#> are not informative for brm_multiple models and are hence not displayed.
#> Please see ?brm_multiple for how to assess convergence in case of such models.
# }
```
