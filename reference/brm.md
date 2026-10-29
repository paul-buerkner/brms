# Fit Bayesian Generalized (Non-)Linear Multivariate Multilevel Models

Fit Bayesian generalized (non-)linear multivariate multilevel models
using Stan for full Bayesian inference. A wide range of distributions
and link functions are supported, allowing users to fit – among others –
linear, robust linear, count data, survival, response times, ordinal,
zero-inflated, hurdle, extended-support beta regression, and even
self-defined mixture models all in a multilevel context. Further
modeling options include non-linear and smooth terms, auto-correlation
structures, censored data, meta-analytic standard errors, and quite a
few more. In addition, all parameters of the response distributions can
be predicted in order to perform distributional regression. Prior
specifications are flexible and explicitly encourage users to apply
prior distributions that actually reflect their beliefs. In addition,
model fit can easily be assessed and compared with posterior predictive
checks and leave-one-out cross-validation.

## Usage

``` r
brm(
  formula,
  data,
  family = gaussian(),
  prior = NULL,
  autocor = NULL,
  data2 = NULL,
  cov_ranef = NULL,
  sample_prior = "no",
  sparse = NULL,
  knots = NULL,
  drop_unused_levels = TRUE,
  stanvars = NULL,
  stan_funs = NULL,
  fit = NA,
  save_pars = getOption("brms.save_pars", NULL),
  save_ranef = NULL,
  save_mevars = NULL,
  save_all_pars = NULL,
  init = NULL,
  inits = NULL,
  chains = 4,
  iter = getOption("brms.iter", 2000),
  warmup = floor(iter/2),
  thin = 1,
  cores = getOption("mc.cores", 1),
  threads = getOption("brms.threads", NULL),
  opencl = getOption("brms.opencl", NULL),
  normalize = getOption("brms.normalize", TRUE),
  control = NULL,
  algorithm = getOption("brms.algorithm", "sampling"),
  backend = getOption("brms.backend", "rstan"),
  future = getOption("future", FALSE),
  silent = getOption("brms.silent", 1),
  seed = NA,
  save_model = NULL,
  stan_model_args = list(),
  file = NULL,
  file_compress = TRUE,
  file_refit = getOption("brms.file_refit", "never"),
  empty = FALSE,
  rename = TRUE,
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

  An object of class `data.frame` (or one that can be coerced to that
  class) containing data of all variables used in the model.

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

- data2:

  A named `list` of objects containing data, which cannot be passed via
  argument `data`. Required for some objects used in autocorrelation
  structures to specify dependency structures as well as for
  within-group covariance matrices.

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

- drop_unused_levels:

  Should unused factors levels in the data be dropped? Defaults to
  `TRUE`.

- stanvars:

  An optional `stanvars` object generated by function
  [`stanvar`](https://paulbuerkner.com/brms/reference/stanvar.md) to
  define additional variables for use in Stan's program blocks.

- stan_funs:

  (Deprecated) An optional character string containing self-defined Stan
  functions, which will be included in the functions block of the
  generated Stan code. It is now recommended to use the `stanvars`
  argument for this purpose instead.

- fit:

  An instance of S3 class `brmsfit` derived from a previous fit;
  defaults to `NA`. If `fit` is of class `brmsfit`, the compiled model
  associated with the fitted result is re-used and all arguments
  modifying the model code or data are ignored. It is not recommended to
  use this argument directly, but to call the
  [`update`](https://paulbuerkner.com/brms/reference/update.brmsfit.md)
  method, instead.

- save_pars:

  An object generated by
  [`save_pars`](https://paulbuerkner.com/brms/reference/save_pars.md)
  controlling which parameters should be saved in the model. The
  argument has no impact on the model fitting itself.

- save_ranef:

  (Deprecated) A flag to indicate if group-level effects for each level
  of the grouping factor(s) should be saved (default is `TRUE`). Set to
  `FALSE` to save memory. The argument has no impact on the model
  fitting itself.

- save_mevars:

  (Deprecated) A flag to indicate if draws of latent noise-free
  variables obtained by using `me` and `mi` terms should be saved
  (default is `FALSE`). Saving these draws allows to better use methods
  such as `predict` with the latent variables but leads to very large R
  objects even for models of moderate size and complexity.

- save_all_pars:

  (Deprecated) A flag to indicate if draws from all variables defined in
  Stan's `parameters` block should be saved (default is `FALSE`). Saving
  these draws is required in order to apply the methods
  `bridge_sampler`, `bayes_factor`, and `post_prob`. Can be set globally
  for the current R session via the `"brms.save_pars"` option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- init:

  Initial values for the sampler. If `NULL` (the default) or `"random"`,
  Stan will randomly generate initial values for parameters in a
  reasonable range. If `0`, all parameters are initialized to zero on
  the unconstrained space. This option is sometimes useful for certain
  families, as it happens that default random initial values cause draws
  to be essentially constant. Generally, setting `init = 0` is worth a
  try, if chains do not initialize or behave well. Alternatively, `init`
  can be a list of lists containing the initial values, or a function
  (or function name) generating initial values. The latter options are
  mainly implemented for internal testing but are available to users if
  necessary. If specifying initial values using a list or a function
  then currently the parameter names must correspond to the names used
  in the generated Stan code (not the names used in R). For more details
  on specifying initial values you can consult the documentation of the
  selected `backend`.

- inits:

  (Deprecated) Alias of `init`.

- chains:

  Number of Markov chains (defaults to 4).

- iter:

  Number of total iterations per chain (including warmup; defaults to
  2000). Can be set globally for the current R session via the
  `"brms.iter"` option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- warmup:

  A positive integer specifying number of warmup (aka burnin)
  iterations. This also specifies the number of iterations used for
  stepsize adaptation, so warmup draws should not be used for inference.
  The number of warmup should not be larger than `iter` and the default
  is `iter/2`.

- thin:

  Thinning rate. Must be a positive integer. Set `thin > 1` to save
  memory and computation time if `iter` is large.

- cores:

  Number of cores to use when executing the chains in parallel, which
  defaults to 1 but we recommend setting the `mc.cores` option to be as
  many processors as the hardware and RAM allow (up to the number of
  chains). For non-Windows OS in non-interactive R sessions, forking is
  used instead of PSOCK clusters.

- threads:

  Number of threads to use in within-chain parallelization. For more
  control over the threading process, `threads` may also be a
  `brmsthreads` object created by
  [`threading`](https://paulbuerkner.com/brms/reference/threading.md).
  Within-chain parallelization is experimental! We recommend its use
  only if you are experienced with Stan's `reduce_sum` function and have
  a slow running model that cannot be sped up by any other means. Can be
  set globally for the current R session via the `"brms.threads"` option
  (see [`options`](https://rdrr.io/r/base/options.html)).

- opencl:

  The platform and device IDs of the OpenCL device to use for fitting
  using GPU support. If you don't know the IDs of your OpenCL device,
  `c(0,0)` is most likely what you need. For more details, see
  [`opencl`](https://paulbuerkner.com/brms/reference/opencl.md). Can be
  set globally for the current R session via the `"brms.opencl"` option

- normalize:

  Logical. Indicates whether normalization constants should be included
  in the Stan code (defaults to `TRUE`). Setting it to `FALSE` requires
  Stan version \>= 2.25 to work. If `FALSE`, sampling efficiency may be
  increased but some post processing functions such as
  [`bridge_sampler`](https://paulbuerkner.com/brms/reference/bridge_sampler.brmsfit.md)
  will not be available. Can be controlled globally for the current R
  session via the \`brms.normalize\` option.

- control:

  A named `list` of parameters to control the sampler's behavior. It
  defaults to `NULL` so all the default values are used. The most
  important control parameters are discussed in the 'Details' section
  below. For a comprehensive overview see
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

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

- backend:

  Character string naming the package to use as the backend for fitting
  the Stan model. Options are `"rstan"` (the default) or `"cmdstanr"`.
  Can be set globally for the current R session via the `"brms.backend"`
  option (see [`options`](https://rdrr.io/r/base/options.html)). Details
  on the rstan and cmdstanr packages are available at
  <https://mc-stan.org/rstan/> and <https://mc-stan.org/cmdstanr/>,
  respectively. Additionally a `"mock"` backend is available to make
  testing brms and packages that depend on it easier. The `"mock"`
  backend does not actually do any fitting, it only checks the generated
  Stan code for correctness and then returns whatever is passed in an
  additional `mock_fit` argument as the result of the fit.

- future:

  Logical; If `TRUE`, the
  [future](https://future.futureverse.org/reference/future.html) package
  is used for parallel execution of the chains and argument `cores` will
  be ignored. Can be set globally for the current R session via the
  `"future"` option. The execution type is controlled via
  [`plan`](https://future.futureverse.org/reference/plan.html) (see the
  examples section below).

- silent:

  Verbosity level between `0` and `2`. If `1` (the default), most of the
  informational messages of compiler and sampler are suppressed. If `2`,
  even more messages are suppressed. The actual sampling progress is
  still printed. Set `refresh = 0` to turn this off as well. If using
  `backend = "rstan"` you can also set `open_progress = FALSE` to
  prevent opening additional progress bars. Can be set globally for the
  current R session via the `"brms.silent"` option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- seed:

  The seed for random number generation to make results reproducible. If
  `NA` (the default), Stan will set the seed randomly.

- save_model:

  Either `NULL` or a character string. In the latter case, the model's
  Stan code is saved via
  [`cat`](https://paulbuerkner.com/brms/reference/addition-terms.md) in
  a text file named after the string supplied in `save_model`.

- stan_model_args:

  A `list` of further arguments passed to
  [`rstan::stan_model`](https://mc-stan.org/rstan/reference/stan_model.html)
  for `backend = "rstan"` or to
  [`cmdstanr::cmdstan_model`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
  for `backend = "cmdstanr"`, which allows to change how models are
  compiled.

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

- empty:

  Logical. If `TRUE`, the Stan model is not created and compiled and the
  corresponding `'fit'` slot of the `brmsfit` object will be empty. This
  is useful if you have estimated a brms-created Stan model outside of
  brms and want to feed it back into the package.

- rename:

  For internal use only.

- ...:

  Further arguments passed to Stan. For `backend = "rstan"` the
  arguments are passed to
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  or
  [`vb`](https://mc-stan.org/rstan/reference/stanmodel-method-vb.html).
  For `backend = "cmdstanr"` the arguments are passed to the
  [`cmdstanr::sample`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)
  or
  [`cmdstanr::variational`](https://mc-stan.org/cmdstanr/reference/model-method-variational.html)
  method.

## Value

An object of class `brmsfit`, which contains the posterior draws along
with many other useful information about the model. Use
`methods(class = "brmsfit")` for an overview on available methods.

## Details

Fit a generalized (non-)linear multivariate multilevel model via full
Bayesian inference using Stan. A general overview is provided in the
vignettes `vignette("brms_overview")` and `vignette("brms_multilevel")`.
For a full list of available vignettes see `vignette(package = "brms")`.

**Formula syntax of brms models**

Details of the formula syntax applied in brms can be found in
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

**Families and link functions**

Details of families supported by brms can be found in
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).

**Prior distributions**

Priors should be specified using the
[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md)
function. Its documentation contains detailed information on how to
correctly specify priors. To find out on which parameters or parameter
classes priors can be defined, use
[`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md).
Default priors are chosen to be non or very weakly informative so that
their influence on the results will be negligible and you usually don't
have to worry about them. However, after getting more familiar with
Bayesian statistics, I recommend you to start thinking about reasonable
informative priors for your model parameters: Nearly always, there is at
least some prior information available that can be used to improve your
inference.

**Adjusting the sampling behavior of Stan**

In addition to choosing the number of iterations, warmup draws, and
chains, users can control the behavior of the NUTS sampler, by using the
`control` argument. The most important reason to use `control` is to
decrease (or eliminate at best) the number of divergent transitions that
cause a bias in the obtained posterior draws. Whenever you see the
warning "There were x divergent transitions after warmup." you should
really think about increasing `adapt_delta`. To do this, write
`control = list(adapt_delta = <x>)`, where `<x>` should usually be value
between `0.8` (current default) and `1`. Increasing `adapt_delta` will
slow down the sampler but will decrease the number of divergent
transitions threatening the validity of your posterior draws.

Another problem arises when the depth of the tree being evaluated in
each iteration is exceeded. This is less common than having divergent
transitions, but may also bias the posterior draws. When it happens,
Stan will throw out a warning suggesting to increase `max_treedepth`,
which can be accomplished by writing
`control = list(max_treedepth = <x>)` with a positive integer `<x>` that
should usually be larger than the current default of `10`. For more
details on the `control` argument see
[`stan`](https://mc-stan.org/rstan/reference/stan.html).

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

## References

Paul-Christian Buerkner (2017). brms: An R Package for Bayesian
Multilevel Models Using Stan. *Journal of Statistical Software*, 80(1),
1-28. `doi:10.18637/jss.v080.i01`

Paul-Christian Buerkner (2018). Advanced Bayesian Multilevel Modeling
with the R Package brms. *The R Journal*. 10(1), 395–411.
`doi:10.32614/RJ-2018-017`

## See also

[`brms`](https://paulbuerkner.com/brms/reference/brms-package.md),
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md),
[`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.md)

## Author

Paul-Christian Buerkner <paul.buerkner@gmail.com>

## Examples

``` r
# \dontrun{
# Poisson regression for the number of seizures in epileptic patients
fit1 <- brm(
  count ~ zBase * Trt + (1|patient),
  data = epilepsy, family = poisson(),
  prior = prior(normal(0, 10), class = b) +
    prior(cauchy(0, 2), class = sd)
)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.32 seconds.
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
#> Chain 1:  Elapsed Time: 1.902 seconds (Warm-up)
#> Chain 1:                1.323 seconds (Sampling)
#> Chain 1:                3.225 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.3e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 2:  Elapsed Time: 1.947 seconds (Warm-up)
#> Chain 2:                1.335 seconds (Sampling)
#> Chain 2:                3.282 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 3:  Elapsed Time: 1.819 seconds (Warm-up)
#> Chain 3:                1.356 seconds (Sampling)
#> Chain 3:                3.175 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.3e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 4:  Elapsed Time: 1.8 seconds (Warm-up)
#> Chain 4:                1.621 seconds (Sampling)
#> Chain 4:                3.421 seconds (Total)
#> Chain 4: 

# generate a summary of the results
summary(fit1)
#>  Family: poisson 
#>   Links: mu = log 
#> Formula: count ~ zBase * Trt + (1 | patient) 
#>    Data: epilepsy (Number of observations: 236) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~patient (Number of levels: 59) 
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(Intercept)     0.59      0.07     0.46     0.73 1.01      718     1517
#> 
#> Regression Coefficients:
#>            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept      1.78      0.12     1.56     2.02 1.01      751     1307
#> zBase          0.70      0.12     0.48     0.93 1.01      712     1058
#> Trt1          -0.29      0.16    -0.62     0.02 1.01      747     1426
#> zBase:Trt1     0.02      0.16    -0.29     0.34 1.01      806     1028
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

# plot the MCMC chains as well as the posterior distributions
plot(fit1)


# predict responses based on the fitted model
head(predict(fit1))
#>      Estimate Est.Error Q2.5 Q97.5
#> [1,]  3.52700  2.087549    0     8
#> [2,]  3.55225  2.084792    0     8
#> [3,]  2.80050  1.897656    0     7
#> [4,]  3.24450  1.929417    0     8
#> [5,] 13.87275  4.141660    6    23
#> [6,]  5.50050  2.626828    1    11

# plot conditional effects for each predictor
plot(conditional_effects(fit1), ask = FALSE)




# investigate model fit
loo(fit1)
#> Warning: Found 8 observations with a pareto_k > 0.7 in model 'fit1'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> 
#> Computed from 4000 by 236 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -670.0 36.2
#> p_loo        92.2 13.6
#> looic      1340.0 72.4
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.5]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     228   96.6%   135     
#>    (0.7, 1]   (bad)        6    2.5%   <NA>    
#>    (1, Inf)   (very bad)   2    0.8%   <NA>    
#> See help('pareto-k-diagnostic') for details.
pp_check(fit1)
#> Using 10 posterior draws for ppc type 'dens_overlay' by default.



# Ordinal regression modeling patient's rating of inhaler instructions
# category specific effects are estimated for variable 'treat'
fit2 <- brm(rating ~ period + carry + cs(treat),
            data = inhaler, family = sratio("logit"),
            prior = set_prior("normal(0,5)"), chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000227 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.27 seconds.
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
#> Chain 1:  Elapsed Time: 2.407 seconds (Warm-up)
#> Chain 1:                2.132 seconds (Sampling)
#> Chain 1:                4.539 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.00021 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.1 seconds.
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
#> Chain 2:  Elapsed Time: 2.474 seconds (Warm-up)
#> Chain 2:                2.121 seconds (Sampling)
#> Chain 2:                4.595 seconds (Total)
#> Chain 2: 
summary(fit2)
#>  Family: sratio 
#>   Links: mu = logit 
#> Formula: rating ~ period + carry + cs(treat) 
#>    Data: inhaler (Number of observations: 572) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept[1]     0.55      0.09     0.38     0.73 1.00     2693     1261
#> Intercept[2]     2.36      0.29     1.86     2.96 1.00     1697     1376
#> Intercept[3]     0.49      0.61    -0.76     1.68 1.00     1664     1290
#> period           0.22      0.16    -0.10     0.54 1.00     2808     1590
#> carry           -0.22      0.16    -0.55     0.09 1.00     1792     1622
#> treat[1]        -0.78      0.23    -1.25    -0.34 1.00     1773     1616
#> treat[2]        -1.02      0.59    -2.24     0.05 1.00     1497     1345
#> treat[3]         1.14      1.22    -1.15     3.56 1.00     1622     1508
#> 
#> Further Distributional Parameters:
#>      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> disc     1.00      0.00     1.00     1.00   NA       NA       NA
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(fit2, ask = FALSE)


WAIC(fit2)
#> Warning: 
#> 4 (0.7%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> 
#> Computed from 2000 by 572 log-likelihood matrix.
#> 
#>           Estimate   SE
#> elpd_waic   -459.5 17.6
#> p_waic         8.2  1.3
#> waic         919.0 35.1
#> 
#> 4 (0.7%) p_waic estimates greater than 0.4. We recommend trying loo instead. 


# Survival regression modeling the time between the first
# and second recurrence of an infection in kidney patients.
fit3 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
            data = kidney, family = lognormal())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: Rejecting initial value:
#> Chain 1:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 1:   Stan can't start sampling from this initial value.
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.25 seconds.
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
#> Chain 1:  Elapsed Time: 0.863 seconds (Warm-up)
#> Chain 1:                0.312 seconds (Sampling)
#> Chain 1:                1.175 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: Rejecting initial value:
#> Chain 2:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.15 seconds.
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
#> Chain 2:  Elapsed Time: 0.912 seconds (Warm-up)
#> Chain 2:                0.345 seconds (Sampling)
#> Chain 2:                1.257 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.6e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
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
#> Chain 3:  Elapsed Time: 0.893 seconds (Warm-up)
#> Chain 3:                0.348 seconds (Sampling)
#> Chain 3:                1.241 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1.6e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
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
#> Chain 4:  Elapsed Time: 0.92 seconds (Warm-up)
#> Chain 4:                0.332 seconds (Sampling)
#> Chain 4:                1.252 seconds (Total)
#> Chain 4: 
summary(fit3)
#>  Family: lognormal 
#>   Links: mu = identity 
#> Formula: time | cens(censored) ~ age * sex + disease + (1 | patient) 
#>    Data: kidney (Number of observations: 76) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~patient (Number of levels: 38) 
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(Intercept)     0.40      0.25     0.02     0.92 1.01     1087     1543
#> 
#> Regression Coefficients:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept         2.62      0.95     0.75     4.45 1.00     2915     2901
#> age               0.02      0.02    -0.03     0.06 1.00     2537     2914
#> sexfemale         2.58      1.13     0.41     4.83 1.00     2627     2589
#> diseaseGN        -0.42      0.51    -1.46     0.60 1.00     3136     2730
#> diseaseAN        -0.53      0.49    -1.49     0.43 1.00     3456     3129
#> diseasePKD        0.59      0.70    -0.80     1.94 1.00     3568     2747
#> age:sexfemale    -0.02      0.03    -0.08     0.02 1.00     2521     2516
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     1.16      0.13     0.93     1.43 1.00     2267     2953
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(fit3, ask = FALSE)


plot(conditional_effects(fit3), ask = FALSE)






# Probit regression using the binomial family
ntrials <- sample(1:10, 100, TRUE)
success <- rbinom(100, size = ntrials, prob = 0.4)
x <- rnorm(100)
data4 <- data.frame(ntrials, success, x)
fit4 <- brm(success | trials(ntrials) ~ x, data = data4,
            family = binomial("probit"))
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.29 seconds.
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
#> Chain 1:  Elapsed Time: 0.073 seconds (Warm-up)
#> Chain 1:                0.072 seconds (Sampling)
#> Chain 1:                0.145 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
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
#> Chain 2:  Elapsed Time: 0.075 seconds (Warm-up)
#> Chain 2:                0.083 seconds (Sampling)
#> Chain 2:                0.158 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.9e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
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
#> Chain 3:  Elapsed Time: 0.075 seconds (Warm-up)
#> Chain 3:                0.069 seconds (Sampling)
#> Chain 3:                0.144 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1.9e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
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
#> Chain 4:  Elapsed Time: 0.074 seconds (Warm-up)
#> Chain 4:                0.081 seconds (Sampling)
#> Chain 4:                0.155 seconds (Total)
#> Chain 4: 
summary(fit4)
#>  Family: binomial 
#>   Links: mu = probit 
#> Formula: success | trials(ntrials) ~ x 
#>    Data: data4 (Number of observations: 100) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept    -0.21      0.05    -0.31    -0.10 1.00     3845     2579
#> x            -0.07      0.05    -0.17     0.02 1.00     3538     2657
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).


# Non-linear Gaussian model
fit5 <- brm(
  bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
     ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
     nl = TRUE),
  data = loss, family = gaussian(),
  prior = c(
    prior(normal(5000, 1000), nlpar = "ult"),
    prior(normal(1, 2), nlpar = "omega"),
    prior(normal(45, 10), nlpar = "theta")
  ),
  control = list(adapt_delta = 0.9)
)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.33 seconds.
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
#> Chain 1:  Elapsed Time: 2.083 seconds (Warm-up)
#> Chain 1:                1.169 seconds (Sampling)
#> Chain 1:                3.252 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
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
#> Chain 2:  Elapsed Time: 2.375 seconds (Warm-up)
#> Chain 2:                1.165 seconds (Sampling)
#> Chain 2:                3.54 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 3:  Elapsed Time: 2.4 seconds (Warm-up)
#> Chain 3:                1.223 seconds (Sampling)
#> Chain 3:                3.623 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.2e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
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
#> Chain 4:  Elapsed Time: 2.344 seconds (Warm-up)
#> Chain 4:                1.079 seconds (Sampling)
#> Chain 4:                3.423 seconds (Total)
#> Chain 4: 
summary(fit5)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: cum ~ ult * (1 - exp(-(dev/theta)^omega)) 
#>          ult ~ 1 + (1 | AY)
#>          omega ~ 1
#>          theta ~ 1
#>    Data: loss (Number of observations: 55) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~AY (Number of levels: 10) 
#>                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(ult_Intercept)   753.65    234.45   424.14  1359.67 1.00     1025     1746
#> 
#> Regression Coefficients:
#>                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> ult_Intercept    5298.39    286.34  4717.95  5873.13 1.01     1006     1387
#> omega_Intercept     1.34      0.05     1.25     1.43 1.00     2308     2611
#> theta_Intercept    46.16      2.05    42.60    50.65 1.00     2202     2160
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma   139.52     15.21   113.51   173.33 1.00     3088     2793
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
conditional_effects(fit5)



# Normal model with heterogeneous variances
data_het <- data.frame(
  y = c(rnorm(50), rnorm(50, 1, 2)),
  x = factor(rep(c("a", "b"), each = 50))
)
fit6 <- brm(bf(y ~ x, sigma ~ 0 + x), data = data_het)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.2 seconds.
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
#> Chain 1:  Elapsed Time: 0.091 seconds (Warm-up)
#> Chain 1:                0.093 seconds (Sampling)
#> Chain 1:                0.184 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.15 seconds.
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
#> Chain 2:  Elapsed Time: 0.097 seconds (Warm-up)
#> Chain 2:                0.091 seconds (Sampling)
#> Chain 2:                0.188 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.5e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.15 seconds.
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
#> Chain 3:  Elapsed Time: 0.098 seconds (Warm-up)
#> Chain 3:                0.085 seconds (Sampling)
#> Chain 3:                0.183 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1.5e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.15 seconds.
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
#> Chain 4:  Elapsed Time: 0.096 seconds (Warm-up)
#> Chain 4:                0.1 seconds (Sampling)
#> Chain 4:                0.196 seconds (Total)
#> Chain 4: 
summary(fit6)
#>  Family: gaussian 
#>   Links: mu = identity; sigma = log 
#> Formula: y ~ x 
#>          sigma ~ 0 + x
#>    Data: data_het (Number of observations: 100) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     0.04      0.13    -0.21     0.28 1.00     5826     3075
#> xb            0.84      0.36     0.16     1.55 1.00     2785     2570
#> sigma_xa     -0.11      0.10    -0.31     0.10 1.00     3243     2313
#> sigma_xb      0.84      0.10     0.66     1.05 1.00     3604     2407
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(fit6)

conditional_effects(fit6)


# extract estimated residual SDs of both groups
sigmas <- exp(as.data.frame(fit6, variable = "^b_sigma_", regex = TRUE))
ggplot(stack(sigmas), aes(values)) +
  geom_density(aes(fill = ind))
#> Error in ggplot(stack(sigmas), aes(values)): could not find function "ggplot"


# Quantile regression predicting the 25%-quantile
fit7 <- brm(bf(y ~ x, quantile = 0.25), data = data_het,
            family = asym_laplace())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.34 seconds.
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
#> Chain 1:  Elapsed Time: 0.196 seconds (Warm-up)
#> Chain 1:                0.209 seconds (Sampling)
#> Chain 1:                0.405 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.32 seconds.
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
#> Chain 2:  Elapsed Time: 0.202 seconds (Warm-up)
#> Chain 2:                0.188 seconds (Sampling)
#> Chain 2:                0.39 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 3.1e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.31 seconds.
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
#> Chain 3:  Elapsed Time: 0.202 seconds (Warm-up)
#> Chain 3:                0.197 seconds (Sampling)
#> Chain 3:                0.399 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3.3e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.33 seconds.
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
#> Chain 4:  Elapsed Time: 0.201 seconds (Warm-up)
#> Chain 4:                0.209 seconds (Sampling)
#> Chain 4:                0.41 seconds (Total)
#> Chain 4: 
summary(fit7)
#>  Family: asym_laplace 
#>   Links: mu = identity 
#> Formula: y ~ x 
#>          quantile = 0.25
#>    Data: data_het (Number of observations: 100) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept    -0.58      0.17    -0.95    -0.26 1.00     3873     2315
#> xb           -0.36      0.32    -0.97     0.28 1.00     2507     2480
#> 
#> Further Distributional Parameters:
#>          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma        0.51      0.05     0.42     0.62 1.00     2555     2532
#> quantile     0.25      0.00     0.25     0.25   NA       NA       NA
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
conditional_effects(fit7)


# use the future package for more flexible parallelization
library(future)
plan(multisession, workers = 4)
fit7a <- update(fit7, future = TRUE)
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 1:  Elapsed Time: 0.203 seconds (Warm-up)
#> Chain 1:                0.191 seconds (Sampling)
#> Chain 1:                0.394 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 2:  Elapsed Time: 0.205 seconds (Warm-up)
#> Chain 2:                0.201 seconds (Sampling)
#> Chain 2:                0.406 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 3.8e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.38 seconds.
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
#> Chain 3:  Elapsed Time: 0.205 seconds (Warm-up)
#> Chain 3:                0.191 seconds (Sampling)
#> Chain 3:                0.396 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3.9e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 4:  Elapsed Time: 0.195 seconds (Warm-up)
#> Chain 4:                0.185 seconds (Sampling)
#> Chain 4:                0.38 seconds (Total)
#> Chain 4: 

# we may also use mirai allowing the use of remote machines as
# well, including the use of computer clusters managed by queuing
# systems
mirai::daemons(4) # distributed computing with url and remote argument
plan(future.mirai::mirai_cluster)
fit7b <- update(fit7, future = TRUE)
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.53 seconds.
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
#> Chain 1:  Elapsed Time: 0.28 seconds (Warm-up)
#> Chain 1:                0.279 seconds (Sampling)
#> Chain 1:                0.559 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.5 seconds.
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
#> Chain 2:  Elapsed Time: 0.261 seconds (Warm-up)
#> Chain 2:                0.189 seconds (Sampling)
#> Chain 2:                0.45 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 5.4e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.54 seconds.
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
#> Chain 3:  Elapsed Time: 0.271 seconds (Warm-up)
#> Chain 3:                0.255 seconds (Sampling)
#> Chain 3:                0.526 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4.1e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.41 seconds.
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
#> Chain 4:  Elapsed Time: 0.218 seconds (Warm-up)
#> Chain 4:                0.18 seconds (Sampling)
#> Chain 4:                0.398 seconds (Total)
#> Chain 4: 
mirai::daemons(0) # shuts down the 4 workers

# fit a model manually via rstan
scode <- stancode(count ~ Trt, data = epilepsy)
sdata <- standata(count ~ Trt, data = epilepsy)
stanfit <- rstan::stan(model_code = scode, data = sdata)
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
#> Chain 1:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 1:                0.012 seconds (Sampling)
#> Chain 1:                0.033 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 2:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 2:                0.011 seconds (Sampling)
#> Chain 2:                0.032 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 3e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 3:  Elapsed Time: 0.022 seconds (Warm-up)
#> Chain 3:                0.012 seconds (Sampling)
#> Chain 3:                0.034 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 4:  Elapsed Time: 0.022 seconds (Warm-up)
#> Chain 4:                0.012 seconds (Sampling)
#> Chain 4:                0.034 seconds (Total)
#> Chain 4: 
# feed the Stan model back into brms
fit8 <- brm(count ~ Trt, data = epilepsy, empty = TRUE)
fit8$fit <- stanfit
fit8 <- rename_pars(fit8)
summary(fit8)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: count ~ Trt 
#>    Data: epilepsy (Number of observations: 236) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     8.46      1.17     6.18    10.69 1.00     3597     2942
#> Trt1         -0.64      1.62    -3.77     2.55 1.00     3568     3139
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma    12.37      0.56    11.28    13.54 1.00     3843     2998
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
# }
```
