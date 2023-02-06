# brms 2.18.0++

### New Features

* Model unstructured autocorrelation matrices via the `unstr` term
thanks to the help of Sebastian Weber. (#1435)
* Model ordinal data with an extra category (non-response or similar)
via the `hurdle_cumulative` family thanks to Stephen Wild. (#1448)
* Improve user control over model recompilation via argument `recompile`
in post-processing methods that require a compiled Stan model.
* Extend control over the `point_estimate` feature in `prepare_predictions`
via the new argument `ndraws_point_estimate`.
* Add support for the latent projection available in **projpred** versions >
2.3.0.

### Bug Fixes

* Fix a Stan syntax error in threaded models with `lasso` priors. (#1427)
* Fix Stan compilation issues for some of the more special 
link functions such as `cauchit` or `softplus`.


# brms 2.18.0

### New Features

* Support regression splines with fixed degrees of freedom 
specified via `s(..., fx = TRUE)`.
* Reuse user-specified control arguments originally passed 
to the Stan backend in `update` and related methods. (#1373, #1378)
* Allow to retain unused factors levels via `drop_unused_levels = FALSE` 
in `brm` and related functions. (#1346)
* Automatically update old default priors based on new input when 
when updating models via `update.brmsfit`. (#1380)
* Allow to use `dirichlet` priors for more parameter types. (#1165)

### Other Changes

* Improve efficiency of converting models fitted with `backend = "cmdstanr"`
to `stanfit` objects thanks to Simon Mills and Jacob Socolar. (#1331)
* Allow for more `O1` optimization of brms-generated Stan models
thanks to Aki Vehtari. (#1382)

### Bug Fixes

* Fix problems with missing boundaries of `sdme` parameters in models
with known response standard errors thanks to Solomon Kurz. (#1348)
* Fix Stan code of `gamma` models with `softplus` link.
* Allow for more flexible data inputs to `brm_multiple`. (#1383)
* Ensure that `control_params` returns the right values for
models fitted with the `cmdstanr` backend. (#1390)
* Fix problems in multivariate spline models when using 
the `subset` addition term. (#1385)


# brms 2.17.0

### New Features

* Add full user control for boundaries of most parameters via the `lb` and 
`ub` arguments of `set_prior` and related functions. (#878, #1094)
* Add family `logistic_normal` for simplex responses. (#1274)
* Add argument `future_args` to `kfold` and `reloo` for additional
control over parallel execution via futures.
* Add families `beta_binomial` & `zero_inflated_beta_binomial` for potentially
over-dispersed and zero-inflated binomial response models thanks to Hayden
Rabel. (#1319 & #1311)
* Display `ppd_*` plots in `pp_check` via argument `prefix`. (#1313)
* Support the `log` link in binomial and beta type families. (#1316)

### Other changes

* Argument `brms_seed` has been added to `get_refmodel.brmsfit()`. (#1287)
* Deprecate argument `inits` in favor of `init` for consistency
with the Stan backends.
* Improve speed of the `summary` method for high-dimensional models. (#1330)

### Bug Fixes

* Fix Stan code of threaded multivariate models
thanks to Anirban Mukherjee. (#1277)
* Fix usage of `int_conditions` in `conditional_smooths`
thanks to Urs Kalbitzer. (#1280)
* Fix an error sometimes occurring for multilevel (reference) models in
`projpred`'s K-fold CV. (#1286)
* Fix response values in `make_standata` for `bernoulli` families
when only 1s are present thanks to Facundo Munoz. (#1298)
* Fix `pp_check` for censored responses to work for all plot types
thanks to Hayden Rabel. (#1327) 
* Ensure that argument `overwrite` in `add_criterion` works as expected 
for all criteria thanks to Andrew Milne. (#1323)
* Fix a problem in `launch_shinystan` occurring when warmup draws 
were saved thanks to Frank Weber. (#1257, #1329)
* Fix numerical stability problems in `log_lik` for ordinal models. (#1192)


# brms 2.16.3

### Other changes

* Move `projpred` from `Imports:` to `Suggests:`. This has the important
implication that users need to load or attach `projpred` themselves if they want
to use it (the more common case is probably attaching, which is achieved by
`library(projpred)`). (#1222)

### Bug Fixes

* Ensure that argument `overwrite` in `add_criterion`
is working as intended thanks to Ruben Arslan. (#1219)
* Fix a bug in `get_refmodel.brmsfit()` (i.e., when using `projpred` for a
`"brmsfit"`) causing offsets not to be recognized. (#1220)
* Several further minor bug fixes.


# brms 2.16.1

### Bug Fixes

* Fix a bug causing problems during post-processing of models
fitted with older versions of brms and the `cmdstanr` backend
thanks to Riccardo Fusaroli. (#1218)


# brms 2.16.0

### New Features

* Support several methods of the `posterior` package. (#1204)
* Substantially extend compatibility of `brms` models
with `emmeans` thanks to Mattan S. Ben-Shachar. (#907, #1134)
* Combine missing value (`mi`) terms with `subset` addition terms. (#1063)
* Expose function `get_dpar` for use in the post-processing
of custom families thank to Martin Modrak. (#1131)
* Support the `squareplus` link function in all families and
distributional parameters that also allow for the `log` link function.
* Add argument `incl_thres` to `posterior_linpred.brmsfit()` allowing to
subtract the threshold-excluding linear predictor from the thresholds in case
of an ordinal family. (#1137)
* Add a `"mock"` backend option to facilitate testing
thanks to Martin Modrak. (#1116)
* Add option `file_refit = "always"` to always overwrite models
stored via the `file` argument. (#1151)
* Initial GPU support via OpenCL thanks to the help
Rok Češnovar. (#1166)
* Support argument `robust` in method `hypothesis`. (#1170)
* Vectorize the Stan code of custom likelihoods via
argument `loop` of `custom_family`. (#1084)
* Experimentally allow category specific effects for
ordinal `cumulative` models. (#1060)
* Regenerate Stan code of an existing model via
argument `regenerate` of method `stancode`.
* Support `expose_functions` for models fitted with the
`cmdstanr` backend thanks to Sebastian Weber. (#1176)
* Support `log_prob` and related functionality in models fitted
with the `cmdstanr` backend via function `add_rstan_model`. (#1184)

### Other Changes

* Remove use of `cbind` to express multivariate models after
over two years of deprecation (please use `mvbind` instead).
* Method `posterior_linpred(transform = TRUE)` is now equal
to `posterior_epred(dpar = "mu")` and no longer deprecated.
* Refactor and extend internal post-processing functions
for ordinal and categorical models thanks to Frank Weber. (#1159)
* Ignore `NA` values in interval censored boundaries as long as
they are unused. (#1070)
* Take offsets into account when deriving default priors
for overall intercept parameters. (#923)
* Soft deprecate measurement error (`me`) terms in favor
of the more general and consistent missing value (`mi`) terms. (#698)

### Bug Fixes

* Fix an issue in the post-processing of non-normal ARMA models
thanks to Thomas Buehrens. (#1149)
* Fix an issue with default baseline hazard knots in `cox` models
thanks to Malcolm Gillies. (#1143)
* Fix a bug in non-linear models caused by accidental
merging of operators in the non-linear formula
thanks to Fernando Miguez. (#1142)
* Correctly trigger a refit for `file_refit = "on_change"` if factor level
names have changed thanks to Martin Modrak. (#1128)
* Validate factors in `validate_newdata` even when they are simultaneously
used as predictors and grouping variables thanks to Martin Modrak. (#1141)
* Fix a bug in the Stan code generation of threaded mixture models
with predicted mixture probabilities thanks to Riccardo Fusaroli. (#1150)
* Remove duplicated Stan code related to the `horseshoe` prior
thanks to Max Joseph. (#1167)
* Fix an issue in the post-processing of non-looped non-linear
parameters thanks to Sebastian Weber.
* Fix an issue in the Stan code of threaded non-looped non-linear
models thanks to Sebastian Weber. (#1175)
* Fix problems in the post-processing of multivariate meta-analytic
models that could lead to incorrect handling of known standard errors.


# brms 2.15.0

### New Features

* Turn off normalization in the Stan model via argument `normalize`.
to increase sampling efficiency thanks to Andrew Johnson. (#1017, #1053)
* Enable `posterior_predict` for truncated continuous models
even if the required CDF or quantile functions are unavailable.
* Update and export `validate_prior` to validate priors supplied by the user.
* Add support for within-chain threading with `rstan (Stan >= 2.25)` backend.
* Apply the R2-D2 shrinkage prior to population-level coefficients
via function `R2D2` to be used in `set_prior`.
* Extend support for `arma` correlation structures in non-normal families.
* Extend scope of variables passed via `data2` for use in the
evaluation of most model terms.
* Refit models previously stored on disc only when necessary thanks to
Martin Modrak. The behavior can be controlled via `file_refit`. (#1058)
* Allow for a finer tuning of informational messages printed in `brm`
via the `silent` argument. (#1076)
* Allow `stanvars` to alter distributional parameters. (#1061)
* Allow `stanvars` to be used inside threaded likelihoods. (#1111)

### Other Changes

* Improve numerical stability of ordinal sequential models
(families `sratio` and `cratio`) thanks to Andrew Johnson. (#1087)

### Bug Fixes

* Allow fitting `multinomial` models with the
`cmdstanr` backend thanks to Andrew Johnson. (#1033)
* Allow user-defined Stan functions in threaded models. (#1034)
* Allow usage of the `:` operator in autocorrelation terms.
* Fix Stan code generation when specifying coefficient-level
priors on spline terms.
* Fix numerical issues occurring in edge cases during
post-processing of Gaussian processes thanks to Marta Kołczyńska.
* Fix an error during post-processing of new levels in
multi-membership terms thanks to Guilherme Mohor.
* Fix a bug in the Stan code of threaded `wiener` drift diffusion
models thanks to the GitHub user yanivabir. (#1085)
* Fix a bug in the threaded Stan code for GPs with categorical
`by` variables thanks to Reece Willoughby. (#1081)
* Fix a bug in the threaded Stan code when using QR decomposition
thanks to Steve Bronder. (#1086)
* Include offsets in `emmeans` related methods thanks to
Russell V. Lenth. (#1096)


# brms 2.14.4

### New Features

* Support `projpred` version 2.0 for variable selection in generalized
linear and additive multilevel models thanks to Alejandro Catalina.
* Support `by` variables in multi-membership terms.
* Use Bayesian bootstrap in `loo_R2`.

### Bug Fixes

* Allow non-linear terms in threaded models.
* Allow multi-membership terms in threaded models.
* Allow `se` addition terms in threaded models.
* Allow `categorical` families in threaded models.
* Fix updating of parameters in `loo_moment_match`.
* Fix facet labels in `conditional_effects` thanks
to Isaac Petersen. (#1014)


# brms 2.14.0

### New Features

* Experimentally support within-chain parallelization via `reduce_sum`
using argument `threads` in `brm` thanks to Sebastian Weber. (#892)
* Add algorithm `fixed_param` to sample from fixed parameter values. (#973)
* No longer remove `NA` values in `data` if there are unused because of
the `subset` addition argument. (#895)
* Combine `by` variables and within-group correlation matrices
in group-level terms. (#674)
* Add argument `robust` to the `summary` method. (#976)
* Parallelize evaluation of the `posterior_predict` and `log_lik`
methods via argument `cores`. (#819)
* Compute effective number of parameters in `kfold`.
* Show prior sources and vectorization in the `print` output
of `brmsprior` objects. (#761)
* Store unused variables in the model's data frame via
argument `unused` of function `brmsformula`.
* Support posterior mean predictions in `emmeans` via
`dpar = "mean"` thanks to Russell V. Lenth. (#993)
* Improve control of which parameters should be saved via
function `save_pars` and corresponding argument in `brm`. (#746)
* Add method `posterior_smooths` to computing predictions
of individual smooth terms. (#738)
* Allow to display grouping variables in `conditional_effects`
using the `effects` argument. (#1012)

### Other Changes

* Improve sampling efficiency for a lot of models by using Stan's
GLM-primitives even in non-GLM cases. (#984)
* Improve sampling efficiency of multilevel models with
within-group covariances thanks to David Westergaard. (#977)
* Deprecate argument `probs` in the `conditional_effects` method
in favor of argument `prob`.

### Bug Fixes

* Fix a problem in `pp_check` inducing wronger observation
orders in time series models thanks to Fiona Seaton. (#1007)
* Fix multiple problems with `loo_moment_match` that prevented
it from working for some more complex models.


# brms 2.13.5

### New Features

* Support the Cox proportional hazards model for
time-to-event data via family `cox`. (#230, #962)
* Support method `loo_moment_match`, which can be used to
update a `loo` object when Pareto k estimates are large.

### Other Changes

* Improve the prediction behavior in post-processing methods
when sampling new levels of grouping factors via
`sample_new_levels = "uncertainty"`. (#956)

### Bug Fixes

* Fix minor problems with MKL on CRAN.


# brms 2.13.3

### New Features

* Fix shape parameters across multiple monotonic terms via argument
`id` in function `mo` to ensure conditionally monotonic effects. (#924)
* Support package `rtdists` as additional backend of `wiener`
distribution functions thanks to the help of Henrik Singmann. (#385)

### Bug Fixes

* Fix generated Stan Code of models with improper global priors and
`constant` priors on some coefficients thanks to Frank Weber. (#919)
* Fix a bug in `conditional_effects` occurring for categorical
models with matrix predictors thanks to Jamie Cranston. (#933)

### Other Changes

* Adjust behavior of the `rate` addition term so that it also
affects the `shape` parameter in `negbinomial` models thanks to
Edward Abraham. (#915)
* Adjust the default inverse-gamma prior on length-scale parameters
of Gaussian processes to be less extreme in edge cases thanks
to Topi Paananen.


# brms 2.13.0

### New Features

* Constrain ordinal thresholds to sum to zero via argument
`threshold` in ordinal family functions thanks to the help of
Marta Kołczyńska.
* Support `posterior_linpred` as method in `conditional_effects`.
* Use `std_normal` in the Stan code for improved efficiency.
* Add arguments `cor`, `id`, and `cov` to the functions `gr` and
`mm` for easy specification of group-level correlation structures.
* Improve workflow to feed back brms-created models which were
fitted somewhere else back into brms. (#745)
* Improve argument `int_conditions` in `conditional_effects` to
work for all predictors not just interactions.
* Support multiple imputation of data passed via `data2` in
`brm_multiple`. (#886)
* Fully support the `emmeans` package thanks to the help
of Russell V. Lenth. (#418)
* Control the within-block position of Stan code added via
`stanvar` using the `position` argument.

### Bug Fixes

* Fix issue in Stan code of models with multiple `me` terms
thanks to Chris Chatham. (#855, #856)
* Fix scaling problems in the estimation of ordinal models with
multiple threshold vectors thanks to Marta Kołczyńska and
Rok Češnovar.
* Allow usage of `std_normal` in `set_prior` thanks to Ben Goodrich. (#867)
* Fix Stan code of distributional models with `weibull`, `frechet`,
or `inverse.gaussian` families thanks to Brian Huey and Jack Caster. (#879)
* Fix Stan code of models which are truncated and weighted at the
same time thanks to Michael Thompson. (#884)
* Fix Stan code of multivariate models with custom families and
data variables passed to the likelihood thanks to Raoul Wolf. (#906)

### Other Changes

* Reduce minimal scale of several default priors from 10 to 2.5.
The resulting priors should remain weakly informative.
* Automatically group observations in `gp` for increased efficiency.
* Rename `parse_bf` to `brmsterms` and deprecate the former function.
* Rename `extract_draws` to `prepare_predictions` and deprecate
the former function.
* Deprecate using a model-dependent `rescor` default.
* Deprecate argument `cov_ranef` in `brm` and related functions.
* Improve several internal interfaces. This should not have any
user-visible changes.
* Simplify the parameterization of the horseshoe prior thanks
to Aki Vehtari. (#873)
* Store fixed distributional parameters as regular draws so that
they behave as if they were estimated in post-processing methods.


# brms 2.12.0

### New Features

* Fix parameters to constants via the `prior` argument. (#783)
* Specify autocorrelation terms directly in the model formula. (#708)
* Translate integer covariates in non-linear formulas to integer
arrays in Stan.
* Estimate `sigma` in combination with fixed correlation matrices
via autocorrelation term `fcor`.
* Use argument `data2` in `brm` and related functions to pass
data objects which cannot be passed via `data`. The usage of `data2`
will be extended in future versions.
* Compute pointwise log-likelihood values via `log_lik` for
non-factorizable Student-t models. (#705)

### Bug Fixes

* Fix output of `posterior_predict` for `multinomial` models
thanks to Ivan Ukhov.
* Fix selection of group-level terms via `re_formula` in
multivariate models thanks to Maxime Dahirel. (#834)
* Enforce correct ordering of terms in `re_formula`
thanks to @ferberkl. (#844)
* Fix post-processing of multivariate multilevel models
when multiple IDs are used for the same grouping factor
thanks to @lott999. (#835)
* Store response category names of ordinal models in the
output of `posterior_predict` again thanks to Mattew Kay. (#838)
* Handle `NA` values more consistently in `posterior_table`
thanks to Anna Hake. (#845)
* Fix a bug in the Stan code of models with multiple monotonic
varying effects across different groups thanks to Julian Quandt.

### Other Changes

* Rename `offset` variables to `offsets` in the generated Stan
code as the former will be reserved in the new stanc3 compiler.


# brms 2.11.1

### Bug Fixes

* Fix version requirement of the `loo` package.
* Fix effective sample size note in the `summary` output. (#824)
* Fix an edge case in the handling of covariates in
special terms thanks to Andrew Milne. (#823)
* Allow restructuring objects multiple times with different
brms versions thanks to Jonathan A. Nations. (#828)
* Fix validation of ordered factors in `newdata`
thanks to Andrew Milne. (#830)

# brms 2.11.0

### New Features

* Support grouped ordinal threshold vectors via addition
argument `resp_thres`. (#675)
* Support method `loo_subsample` for performing approximate
leave-one-out cross-validation for large data.
* Allow storing more model fit criteria via `add_criterion`. (#793)

### Bug Fixes

* Fix prediction uncertainties of new group levels for
`sample_new_levels = "uncertainty"` thanks to Dominic Magirr. (#779)
* Fix problems when using `pp_check` on
censored models thanks to Andrew Milne. (#744)
* Fix error in the generated Stan code of multivariate
`zero_inflated_binomial` models thanks to Raoul Wolf. (#756)
* Fix predictions of spline models when using addition
argument `subset` thanks to Ruben Arslan.
* Fix out-of-sample predictions of AR models when predicting
more than one step ahead.
* Fix problems when using `reloo` or `kfold` with CAR models.
* Fix problems when using `fitted(..., scale = "linear")` with
multinomial models thanks to Santiago Olivella. (#770)
* Fix problems in the `as.mcmc` method for thinned models
thanks to @hoxo-m. (#811)
* Fix problems in parsing covariates of special effects terms
thanks to Riccardo Fusaroli (#813)

### Other Changes

* Rename `marginal_effects` to `conditional_effects` and
`marginal_smooths` to `conditional_smooths`. (#735)
* Rename `stanplot` to `mcmc_plot`.
* Add method `pp_expect` as an alias of `fitted`. (#644)
* Model fit criteria computed via `add_criterion` are now
stored in the `brmsfit$criteria` slot.
* Deprecate `resp_cat` in favor of `resp_thres`.
* Deprecate specifying global priors on regression coefficients
in categorical and multivariate models.
* Improve names of weighting methods in `model_weights`.
* Deprecate reserved variable `intercept` in favor of `Intercept`.
* Deprecate argument `exact_match` in favor of `fixed`.
* Deprecate functions `add_loo` and `add_waic`
in favor of `add_criterion`.

# brms 2.10.0

### New Features

* Improve convergence diagnostics in the `summary` output. (#712)
* Use primitive Stan GLM functions whenever possible. (#703)
* Pass real and integer data vectors to custom families via
the addition arguments `vreal` and `vint`. (#707)
* Model compound symmetry correlations via `cor_cosy`. (#403)
* Predict `sigma` in combination with several
autocorrelation structures. (#403)
* Use addition term `rate` to conveniently handle
denominators of rate responses in log-linear models.
* Fit BYM2 CAR models via `cor_car` thanks to the case study
and help of Mitzi Morris.

### Other Changes

* Substantially improve the sampling efficiency of SAR models
thanks to the GitHub user aslez. (#680)
* No longer allow changing the boundaries
of autocorrelation parameters.
* Set the number of trials to 1 by default in
`marginal_effects` if not specified otherwise. (#718)
* Use non-standard evaluation for addition terms.
* Name temporary intercept parameters more consistently
in the Stan code.

### Bug Fixes

* Fix problems in the post-processing of `me` terms with
grouping factors thanks to the GitHub user tatters. (#706)
* Allow grouping variables to start with a dot
thanks to Bruno Nicenboim. (#679)
* Allow the `horseshoe` prior in categorical and
related models thanks to the Github user tatters. (#678)
* Fix extraction of prior samples for overall intercepts in
`prior_samples` thanks to Jonas Kristoffer Lindelov. (#696)
* Allow underscores to be used in category names
of categorical responses thanks to Emmanuel Charpentier. (#672)
* Fix Stan code of multivariate models with multi-membership
terms thanks to the Stan discourse user Pia.
* Improve checks for non-standard variable names
thanks to Ryan Holbrook. (#721)
* Fix problems when plotting facetted spaghetti plots
via `marginal_smooths` thanks to Gavin Simpson. (#740)


# brms 2.9.0

### New Features

* Specify non-linear ordinal models. (#623)
* Allow to fix thresholds in ordinal mixture models (#626)
* Use the `softplus` link function in various families. (#622)
* Use QR decomposition of design matrices via argument
`decomp` of `brmsformula` thanks to the help of Ben Goodrich. (#640)
* Define argument `sparse` separately for each model formula.
* Allow using `bayes_R2` and `loo_R2` with ordinal models. (#639)
* Support `cor_arma` in non-normal models. (#648)

### Other Changes

* Change the parameterization of monotonic effects to
improve their interpretability. (#578)
* No longer support the `cor_arr` and `cor_bsts` correlation
structures after a year of deprecation.
* Refactor internal evaluation of special predictor terms.
* Improve penalty of splines thanks to Ben Goodrich
and Ruben Arslan.

### Bug Fixes

* Fix a problem when applying `marginal_effects` to
measurement error models thanks to Jonathan A. Nations. (#636)
* Fix computation of log-likelihood values for weighted
mixture models.
* Fix computation of fitted values for truncated lognormal
and weibull models.
* Fix checking of response boundaries for models with
missing values thanks to Lucas Deschamps.
* Fix Stan code of multivariate models with both residual
correlations and missing value terms thanks to Solomon Kurz.
* Fix problems with interactions of special terms
when extracting variable names in `marginal_effects`.
* Allow compiling a model in `brm_multiple` without
sampling thanks to Will Petry. (#671)


# brms 2.8.0

### New Features

* Fit multinomial models via family `multinomial`. (#463)
* Fit Dirichlet models via family `dirichlet`. (#463)
* Fit conditional logistic models using the `categorical` and
`multinomial` families together with non-linear formula syntax. (#560)
* Choose the reference category of `categorical` and related
families via argument `refcat` of the corresponding family functions.
* Use different subsets of the data in different univariate parts
of a multivariate model via addition argument `subset`. (#360)
* Control the centering of population-level design matrices
via argument `center` of `brmsformula` and related functions.
* Add an `update` method for `brmsfit_multiple` objects. (#615)
* Split folds after `group` in the `kfold` method. (#619)

### Other changes

* Deprecate `compare_ic` and instead recommend `loo_compare` for the
comparison of `loo` objects to ensure consistency between packages. (#414)
* Use the **glue** package in the Stan code generation. (#549)
* Introduce `mvbind` to eventually replace `cbind`
in the formula syntax of multivariate models.
* Validate several sampling-related arguments in
`brm` before compiling the Stan model. (#576)
* Show evaluated vignettes on CRAN again. (#591)
* Export function `get_y` which is used to extract response
values from `brmsfit` objects.

### Bug fixes

* Fix an error when trying to change argument `re_formula`
in `bayes_R2` thanks to the GitHub user emieldl. (#592)
* Fix occasional problems when running chains in parallel
via the **future** package thanks to Jared Knowles. (#579)
* Ensure correct ordering of response categories in ordinal
models thanks to Jonas Kristoffer Lindelov. (#580)
* Ignore argument `resp` of `marginal_effects` in
univariate models thanks to Vassilis Kehayas. (#589)
* Correctly disable cell-mean coding in varying effects.
* Allow to fix parameter `ndt` in drift diffusion models.
* Fix Stan code for t-distributed varying effects
thanks to Ozgur Asar.
* Fix an error in the post-processing of monotonic effects
occurring for multivariate models thanks to James Rae. (#598)
* Fix lower bounds in truncated discrete models.
* Fix checks of the original data in `kfold` thanks to
the GitHub user gcolitti. (#602)
* Fix an error when applying the `VarCorr` method to
meta-analytic models thanks to Michael Scharkow. (#616)


# brms 2.7.0

### New features

* Fit approximate and non-isotropic Gaussian processes via `gp`. (#540)
* Enable parallelization of model fitting in `brm_multiple`
via the future package. (#364)
* Perform posterior predictions based on k-fold cross-validation
via `kfold_predict`. (#468)
* Indicate observations for out-of-sample predictions in
ARMA models via argument `oos` of `extract_draws`. (#539)

### Other changes

* Allow factor-like variables in smooth terms. (#562)
* Make plotting of `marginal_effects` more robust to
the usage of non-standard variable names.
* Deactivate certain data validity checks when using custom families.
* Improve efficiency of adjacent category models.
* No longer print informational messages from the Stan parser.

### Bug fixes

* Fix an issue that could result in a substantial efficiency
drop of various post-processing methods for larger models.
* Fix an issue when that resulted in an error when
using `fitted(..., scale = "linear")` with ordinal models
thanks to Andrew Milne. (#557)
* Allow setting priors on the overall intercept in sparse models.
* Allow sampling from models with only a single observation
that also contain an offset thanks to Antonio Vargas. (#545)
* Fix an error when sampling from priors in mixture models
thanks to Jacki Buros Novik. (#542)
* Fix a problem when trying to sample from priors of
parameter transformations.
* Allow using `marginal_smooths` with ordinal models
thanks to Andrew Milne. (#570)
* Fix an error in the post-processing of `me`
terms thanks to the GitHub user hlluik. (#571)
* Correctly update `warmup` samples when using
`update.brmsfit`.


# brms 2.6.0

### New features

* Fit factor smooth interactions thanks to Simon Wood.
* Specify separate priors for thresholds in ordinal models. (#524)
* Pass additional arguments to `rstan::stan_model` via argument
`stan_model_args` in `brm`. (#525)
* Save model objects via argument `file` in `add_ic`
after adding model fit criteria. (#478)
* Compute density ratios based on MCMC samples via `density_ratio`.
* Ignore offsets in various post-processing methods via
argument `offset`.
* Update addition terms in formulas via `update_adterms`.

### Other changes

* Improve internal modularization of smooth terms.
* Reduce size of internal example models.

### Bug fixes

* Correctly plot splines with factorial covariates via `marginal_smooths`.
* Allow sampling from priors in intercept only models
thanks to Emmanuel Charpentier. (#529)
* Allow logical operators in non-linear formulas.


# brms 2.5.0

### New features

* Improve `marginal_effects` to better display ordinal and
categorical models via argument `categorical`. (#491, #497)
* Improve method `kfold` to offer more options for specifying
omitted subsets. (#510)
* Compute estimated values of non-linear parameters via
argument `nlpar` in method `fitted`.
* Disable automatic cell-mean coding in model formulas without
an intercept via argument `cmc` of `brmsformula` and related
functions thanks to Marie Beisemann.
* Allow using the `bridge_sampler` method even if
prior samples are drawn within the model. (#485)
* Specify post-processing functions of custom families
directly in `custom_family`.
* Select a subset of coefficients in `fixef`, `ranef`,
and `coef` via argument `pars`. (#520)
* Allow to `overwrite` already stored fit indices
when using `add_ic`.

### Other changes

* Ignore argument `resp` when post-processing
univariate models thanks to Ruben Arslan. (#488)
* Deprecate argument `ordinal` of `marginal_effects`. (#491)
* Deprecate argument `exact_loo` of `kfold`. (#510)
* Deprecate usage of `binomial` families without specifying `trials`.
* No longer sample from priors of population-level intercepts
when using the default intercept parameterization.

### Bug fixes

* Correctly sample from LKJ correlation priors
thanks to Donald Williams.
* Remove stored fit indices when calling `update` on
brmsfit objects thanks to Emmanuel Charpentier. (#490)
* Fix problems when predicting a single data point using
spline models thanks to Emmanuel Charpentier. (#494)
* Set `Post.Prob = 1` if `Evid.Ratio = Inf` in
method `hypothesis` thanks to Andrew Milne. (#509)
* Ensure correct handling of argument `file` in `brm_multiple`.


# brms 2.4.0

### New features

* Define custom variables in all of Stan's program blocks
via function `stanvar`. (#459)
* Change the scope of non-linear parameters to be global
within univariate models. (#390)
* Allow to automatically group predictor values in Gaussian
processes specified via `gp`. This may lead to a
considerable increase in sampling efficiency. (#300)
* Compute LOO-adjusted R-squared using method `loo_R2`.
* Compute non-linear predictors outside of a loop over
observations by means of argument `loop` in `brmsformula`.
* Fit non-linear mixture models. (#456)
* Fit censored or truncated mixture models. (#469)
* Allow `horseshoe` and `lasso` priors to be set on special
population-level effects.
* Allow vectors of length greater one to be passed to `set_prior`.
* Conveniently save and load fitted model objects in `brm`
via argument `file`. (#472)
* Display posterior probabilities in the output of `hypothesis`.

### Other changes

* Deprecate argument `stan_funs` in `brm` in favor of using the
`stanvars` argument for the specification of custom Stan functions.
* Deprecate arguments `flist` and `...` in `nlf`.
* Deprecate argument `dpar` in `lf` and `nlf`.

### Bug fixes

* Allow custom families in mixture models thanks to Noam Ross. (#453)
* Ensure compatibility with **mice** version 3.0. (#455)
* Fix naming of correlation parameters of group-level terms
with multiple subgroups thanks to Kristoffer Magnusson. (#457)
* Improve scaling of default priors in `lognormal` models (#460).
* Fix multiple problems in the post-processing of categorical models.
* Fix validation of nested grouping factors in post-processing
methods when passing new data thanks to Liam Kendall.


# brms 2.3.1

### New features

* Allow censoring and truncation in zero-inflated and hurdle models. (#430)
* Export zero-inflated and hurdle distribution functions.

### Other changes

* Improve sampling efficiency of the ordinal families
`cumulative`, `sratio`, and `cratio`. (#433)
* Allow to specify a single k-fold subset in method `kfold`. (#441)

### Bug fixes

* Fix a problem in `launch_shinystan` due to which the
maximum treedepth was not correctly displayed thanks to
Paul Galpern. (#431)


# brms 2.3.0

### Features

* Extend `cor_car` to support intrinsic CAR models in pairwise
difference formulation thanks to the case study of Mitzi Morris.
* Compute `loo` and related methods for non-factorizable normal models.

### Other changes

* Rename quantile columns in `posterior_summary`. This affects the
output of `predict` and related methods if `summary = TRUE`. (#425)
* Use hashes to check if models have the same response values
when performing model comparisons. (#414)
* No longer set `pointwise` dynamically in `loo` and related methods. (#416)
* No longer show information criteria in the summary output.
* Simplify internal workflow to implement native response distributions. (#421)

### Bug fixes

* Allow `cor_car` in multivariate models with residual correlations
thanks to Quentin Read. (#427)
* Fix a problem in the Stan code generation of distributional `beta` models
thanks to Hans van Calster. (#404)
* Fix `launch_shinystan.brmsfit` so that all parameters
are now shown correctly in the diagnose tab. (#340)


# brms 2.2.0

### Features

* Specify custom response distributions with function `custom_family`. (#381)
* Model missing values and measurement error in responses using the `mi`
addition term. (#27, #343)
* Allow missing values in predictors using `mi` terms on the right-hand side of
model formulas. (#27)
* Model interactions between the special predictor terms `mo`, `me`, and `mi`.
(#313)
* Introduce methods `model_weights` and `loo_model_weights` providing several
options to compute model weights. (#268)
* Introduce method `posterior_average` to extract posterior samples averaged
across models. (#386)
* Allow hyperparameters of group-level effects to vary over the levels of a
categorical covariate using argument `by` in function `gr`. (#365)
* Allow predictions of measurement-error models with new data. (#335)
* Pass user-defined variables to Stan via `stanvar`. (#219, #357)
* Allow ordinal families in mixture models. (#389)
* Model covariates in multi-membership structures that vary over the levels of
the grouping factor via `mmc` terms. (#353)
* Fit shifted log-normal models via family `shifted_lognormal`. (#218)
* Specify nested non-linear formulas.
* Introduce function `make_conditions` to ease preparation of conditions for
`marginal_effects`.

### Other changes

* Change the parameterization of `weibull` and `exgaussian` models to be
consistent with other model classes. Post-processing of related models fitted
with earlier version of `brms` is no longer possible.
* Treat integer responses in `ordinal` models as directly indicating categories
even if the lowest integer is not one.
* Improve output of the `hypothesis` method thanks to the ideas of Matti Vuorre.
(#362)
* Always plot `by` variables as facets in `marginal_smooths`.
* Deprecate the `cor_bsts` correlation structure.

### Bug fixes

* Allow the `:` operator to combine groups in multi-membership terms thanks to
Gang Chen.
* Avoid an unexpected error when calling `LOO` with argument `reloo = TRUE`
thanks to Peter Konings. (#348)
* Fix problems in `predict` when applied to categorical models thanks to Lydia
Andreyevna Krasilnikova and Thomas Vladeck. (#336, #345)
* Allow truncation in multivariate models with missing values thanks to Malte
Lau Petersen. (#380)
* Force time points to be unique within groups in autocorrelation structures
thanks to Ruben Arslan. (#363)
* Fix problems when post-processing multiple uncorrelated group-level terms of
the same grouping factor thanks to Ivy Jansen. (#374)
* Fix a problem in the Stan code of multivariate `weibull` and `frechet` models
thanks to the GitHub user philj1s. (#375)
* Fix a rare error when post-processing `binomial` models thanks to the GitHub
user SeanH94. (#382)
* Keep attributes of variables when preparing the `model.frame` thanks to Daniel
Luedecke. (#393)


# brms 2.1.0

### Features

* Fit models on multiple imputed datasets via `brm_multiple` thanks to Ruben
Arslan. (#27)
* Combine multiple `brmsfit` objects via function `combine_models`.
* Compute model averaged posterior predictions with method `pp_average`. (#319)
* Add new argument `ordinal` to `marginal_effects` to generate special plots for
ordinal models thanks to the idea  of the GitHub user silberzwiebel. (#190)
* Use informative inverse-gamma priors for length-scale parameters of Gaussian
processes. (#275)
* Compute hypotheses for all levels of a grouping factor at once using
argument `scope` in method `hypothesis`. (#327)
* Vectorize user-defined `Stan` functions exported via
`export_functions` using argument `vectorize`.
* Allow predicting new data in models with ARMA autocorrelation structures.


### Bug fixes

* Correctly recover noise-free coefficients through `me` terms thanks to Ruben
Arslan. As a side effect, it is no longer possible to define priors on
noise-free `Xme` variables directly, but only on their hyper-parameters `meanme`
and `sdme`.
* Fix problems in renaming parameters of the `cor_bsts` structure thanks to
Joshua Edward Morten. (#312)
* Fix some unexpected errors when predicting from ordinal models thanks to David
Hervas and Florian Bader. (#306, #307, #331)
* Fix problems when estimating and predicting multivariate ordinal models thanks
to David West. (#314)
* Fix various minor problems in autocorrelation structures thanks to David West.
(#320)




# brms 2.0.1

### Features

* Export the helper functions `posterior_summary` and `posterior_table` both
being used to summarize posterior samples and predictions.


### Bug fixes

* Fix incorrect computation of intercepts in `acat` and `cratio` models thanks
to Peter Phalen. (#302)
* Fix `pointwise` computation of `LOO` and `WAIC` in multivariate models with
estimated residual correlation structure.
* Fix problems in various S3 methods sometimes requiring unused variables to be
specified in `newdata`.
* Fix naming of Stan models thanks to Hao Ran Lai.




# brms 2.0.0

This is the second major release of `brms`. The main new feature are generalized
multivariate models, which now support everything already possible in univariate
models, but with multiple response variables. Further, the internal structure of
the package has been improved considerably to be easier to maintain and extend
in the future. In addition, most deprecated functionality and arguments have
been removed to provide a clean new start for the package. Models fitted with
`brms` 1.0 or higher should remain fully compatible with `brms` 2.0.

### Features

* Add support for generalized multivariate models, where each of the univariate
models may have a different family and autocorrelation structure. Residual
correlations can be estimated for multivariate `gaussian` and `student` models.
All features supported in univariate models are now also available in
multivariate models. (#3)
* Specify different formulas for different categories in `categorical` models.
* Add weakly informative default priors for the parameter class `Intercept` to
improve convergence of more complex distributional models.
* Optionally display the MC standard error in the `summary` output. (#280)
* Add argument `re.form` as an alias of `re_formula` to the methods
`posterior_predict`, `posterior_linpred`, and `predictive_error` for consistency
with other packages making use of these methods. (#283)


### Other changes

* Refactor many parts of the package to make it more consistent and easier to
extend.
* Show the link functions of all distributional parameters in the `summary`
output. (#277)
* Reduce working memory requirements when extracting posterior samples for use
in `predict` and related methods thanks to Fanyi Zhang. (#224)
* Remove deprecated aliases of functions and arguments from the package. (#278)
* No longer support certain prior specifications, which were previously labeled
as deprecated.
* Remove the deprecated addition term `disp` from the package.
* Remove old versions of methods `fixef`, `ranef`, `coef`, and `VarCorr`.
* No longer support models fitted with `brms` < 1.0, which used the multivariate
`'trait'` syntax originally deprecated in `brms` 1.0.
* Make posterior sample extraction in the `summary` method cleaner and less
error prone.
* No longer fix the seed for random number generation in `brm` to avoid
unexpected behavior in simulation studies.


### Bug fixes

* Store `stan_funs` in `brmsfit` objects to allow using `update` on models with
user-defined Stan functions thanks to Tom Wallis. (#288)
* Fix problems in various post-processing methods when applied to models with
the reserved variable `intercept` in group-level terms thanks to the GitHub user
ASKurz. (#279)
* Fix an unexpected error in `predict` and related methods when setting
`sample_new_levels = "gaussian"` in models with only one group-level effect.
Thanks to Timothy Mastny. (#286)




# brms 1.10.2

### Features

* Allow setting priors on noise-free variables specified via function `me`.
* Add arguments `Ksub`, `exact_loo` and `group` to method `kfold` for defining
omitted subsets according to a grouping variable or factor.
* Allow addition argument `se` in `skew_normal` models.


### Bug fixes

* Ensure correct behavior of horseshoe and lasso priors in multivariate models
thanks to Donald Williams.
* Allow using `identity` links on all parameters of the `wiener` family thanks
to Henrik Singmann. (#276)
* Use reasonable dimnames in the output of `fitted` when returning linear
predictors of ordinal models thanks to the GitHub user atrolle. (#274)
* Fix problems in `marginal_smooths` occurring for multi-membership models
thanks to Hans Tierens.




# brms 1.10.0

### Features

* Rebuild monotonic effects from scratch to allow specifying interactions with
other variables. (#239)
* Introduce methods `posterior_linpred` and `posterior_interval` for consistency
with other model fitting packages based on `Stan`.
* Introduce function `theme_black` providing a black `ggplot2` theme.
* Specify special group-level effects within the same terms as ordinary
group-level effects.
* Add argument `prob` to `summary`, which allows to control the width of the
computed uncertainty intervals. (#259)
* Add argument `newdata` to the `kfold` method.
* Add several arguments to the `plot` method of `marginal_effects` to improve
control over the appearences of the plots.


### Other changes

* Use the same noise-free variables for all model parts in measurement error
models. (#257)
* Make names of local-level terms used in the `cor_bsts` structure more
informative.
* Store the `autocor` argument within `brmsformula` objects.
* Store posterior and prior samples in separate slots in the output of method
`hypothesis`.
* No longer change the default theme of `ggplot2` when attaching `brms`. (#256)
* Make sure signs of estimates are not dropped when rounding to zero in
`summary.brmsfit`. (#263)
* Refactor parts of `extract_draws` and `linear_predictor` to be more consistent
with the rest of the package.


### Bug fixes

* Do not silence the `Stan` parser when calling `brm` to get informative error
messages about invalid priors.
* Fix problems with spaces in priors passed to `set_prior`.
* Handle non `data.frame` objects correctly in `hypothesis.default`.
* Fix a problem relating to the colour of points displayed in
`marginal_effects`.



# brms 1.9.0

### Features

* Perform model comparisons based on marginal likelihoods using the methods
`bridge_sampler`, `bayes_factor`, and `post_prob` all powered by the
`bridgesampling` package.
* Compute a Bayesian version of R-squared with the `bayes_R2` method.
* Specify non-linear models for all distributional parameters.
* Combine multiple model formulas using the `+` operator and the helper
functions `lf`, `nlf`, and `set_nl`.
* Combine multiple priors using the `+` operator.
* Split the `nlpar` argument of `set_prior` into the three arguments `resp`,
`dpar`, and `nlpar` to allow for more flexible prior specifications.


### Other changes

* Refactor parts of the package to prepare for the implementation of more
flexible multivariate models in future updates.
* Keep all constants in the log-posterior in order for `bridge_sampler` to be
working correctly.
* Reduce the amount of renaming done within the `stanfit` object.
* Rename argument `auxpar` of `fitted.brmsfit` to `dpar`.
* Use the `launch_shinystan` generic provided by the `shinystan` package.
* Set `bayesplot::theme_default()` as the default `ggplot2` theme when attaching
`brms`.
* Include citations of the `brms` overview paper as published in the Journal of
Statistical Software.


### Bug fixes

* Fix problems when calling `fitted` with `hurdle_lognormal` models thanks to
Meghna Krishnadas.
* Fix problems when predicting `sigma` in `asym_laplace` models thanks to Anna
Josefine Sorensen.



# brms 1.8.0

### Features

* Fit conditional autoregressive (CAR) models via function `cor_car` thanks to
the case study of Max Joseph.
* Fit spatial autoregressive (SAR) models via function `cor_sar`. Currently
works for families `gaussian` and `student`.
* Implement skew normal models via family `skew_normal`. Thanks to Stephen
Martin for suggestions on the parameterization.
* Add method `reloo` to perform exact cross-validation for problematic
observations and `kfold` to perform k-fold cross-validation thanks to the Stan
Team.
* Regularize non-zero coefficients in the `horseshoe` prior thanks to Juho
Piironen and Aki Vehtari.
* Add argument `new_objects` to various post-processing methods to allow for
passing of data objects, which cannot be passed via `newdata`.
* Improve parallel execution flexibility via the `future` package.


### Other changes

* Improve efficiency and stability of ARMA models.
* Throw an error when the intercept is removed in an ordinal model instead of
silently adding it back again.
* Deprecate argument `threshold` in `brm` and instead recommend passing
`threshold` directly to the ordinal family functions.
* Throw an error instead of a message when invalid priors are passed.
* Change the default value of the `autocor` slot in `brmsfit` objects to an
empty `cor_brms` object.
* Shorten `Stan` code by combining declarations and definitions where possible.


### Bug fixes

* Fix problems in `pp_check` when the variable specified in argument `x` has
attributes thanks to Paul Galpern.
* Fix problems when computing fitted values for truncated discrete models based
on new data thanks to Nathan Doogan.
* Fix unexpected errors when passing models, which did not properly initialize,
to various post-processing methods.
* Do not accidently drop the second dimension of matrices in `summary.brmsfit`
for models with only a single observation.




# brms 1.7.0

### Features

* Fit latent Gaussian processes of one or more covariates via function `gp`
specified in the model formula (#221).
* Rework methods `fixef`, `ranef`, `coef`, and `VarCorr` to be more flexible and
consistent with other post-processing methods (#200).
* Generalize method `hypothesis` to be applicable on all objects coercible to a
`data.frame` (#198).
* Visualize predictions via spaghetti plots using argument `spaghetti` in
`marginal_effects` and `marginal_smooths`.
* Introduce method `add_ic` to store and reuse information criteria in fitted
model objects (#220).
* Allow for negative weights in multi-membership grouping structures.
* Introduce an `as.array` method for `brmsfit` objects.


### Other changes

* Show output of \R code in HTML vignettes thanks to Ben Goodrich (#158).
* Resolve citations in PDF vignettes thanks to Thomas Kluth (#223).
* Improve sampling efficiency for `exgaussian` models thanks to Alex Forrence
(#222).
* Also transform data points when using argument `transform` in
`marginal_effects` thanks to Markus Gesmann.


### Bug fixes

* Fix an unexpected error in `marginal_effects` occurring for some models with
autocorrelation terms thanks to Markus Gesmann.
* Fix multiple problems occurring for models with the `cor_bsts` structure
thanks to Andrew Ellis.



# brms 1.6.1

### Features

* Implement zero-one-inflated beta models via family `zero_one_inflated_beta`.
* Allow for more link functions in zero-inflated and hurdle models.


### Other changes

* Ensure full compatibility with `bayesplot` version 1.2.0.
* Deprecate addition argument `disp`.


### Bug fixes

* Fix problems when setting priors on coefficients of auxiliary parameters when
also setting priors on the corresponding coefficients of the mean parameter.
Thanks to Matti Vuorre for reporting this bug.
* Allow ordered factors to be used as grouping variables thanks to the GitHub
user itissid.



# brms 1.6.0

### Features

* Fit finite mixture models using family function `mixture`.
* Introduce method `pp_mixture` to compute posterior probabilities of mixture
component memberships thanks to a discussion with Stephen Martin.
* Implement different ways to sample new levels of grouping factors in `predict`
and related methods through argument `sample_new_levels`. Thanks to Tom Wallis
and Jonah Gabry for a detailed discussion about this feature.
* Add methods `loo_predict`, `loo_linpred`, and `loo_predictive_interval` for
computing LOO predictions thanks to Aki Vehtari and Jonah Gabry.
* Allow using `offset` in formulas of non-linear and auxiliary parameters.
* Allow sparse matrix multiplication in non-linear and distributional models.
* Allow using the `identity` link for all auxiliary parameters.
* Introduce argument `negative_rt` in `predict` and `posterior_predict` to
distinguish responses on the upper and lower boundary in `wiener` diffusion
models thanks to Guido Biele.
* Introduce method `control_params` to conveniently extract control parameters
of the NUTS sampler.
* Introduce argument `int_conditions` in `marginal_effects` for enhanced
plotting of two-way interactions thanks to a discussion with Thomas Kluth.
* Improve flexibility of the `conditions` argument of `marginal_effects`.
* Extend method `stanplot` to correctly handle some new `mcmc_` plots of the
`bayesplot` package.


### Other changes

* Improve the `update` method to only recompile models when the `Stan` code
changes.
* Warn about divergent transitions when calling `summary` or `print` on
`brmsfit` objects.
* Warn about unused variables in argument `conditions` when calling
`marginal_effects`.
* Export and document several distribution functions that were previously kept
internal.


### Bug fixes

* Fix problems with the inclusion of offsets occurring for more complicated
formulas thanks to Christian Stock.
* Fix a bug that led to invalid Stan code when sampling from priors in intercept
only models thanks to Tom Wallis.
* Correctly check for category specific group-level effects in non-ordinal
models thanks to Wayne Folta.
* Fix problems in `pp_check` when specifying argument `newdata` together with
arguments `x` or `group`.
* Rename the last column in the output of `hypothesis` to `"star"` in order to
avoid problems with zero length column names thanks to the GitHub user
puterleat.
* Add a missing new line statement at the end of the `summary` output thanks to
Thomas Kluth.



# brms 1.5.1

### Features

* Allow `horseshoe` and `lasso` priors to be applied on population-level effects
of non-linear and auxiliary parameters.
* Force recompiling `Stan` models in `update.brmsfit` via argument `recompile`.


### Other changes

* Avoid indexing of matrices in non-linear models to slightly improve sampling
speed.


### Bug fixes

* Fix a severe problem (introduced in version 1.5.0), when predicting `Beta`
models thanks to Vivian Lam.
* Fix problems when summarizing some models fitted with older version of `brms`
thanks to Vivian Lam.
* Fix checks of argument `group` in method `pp_check` thanks to Thomas K.
* Get arguments `subset` and `nsamples` working correctly in `marginal_smooths`.



# brms 1.5.0

### Features

* Implement the generalized extreme value distribution via family
`gen_extreme_value`.
* Improve flexibility of the `horseshoe` prior thanks to Juho Piironen.
* Introduce auxiliary parameter `mu` as an alternative to specifying effects
within the `formula` argument in function `brmsformula`.
* Return fitted values of auxiliary parameters via argument `auxpar` of method
`fitted`.
* Add vignette `"brms_multilevel"`, in which the advanced formula syntax of
`brms` is explained in detail using several examples.


### Other changes

* Refactor various parts of the package to ease implementation of mixture and
multivariate models in future updates. This should not have any user visible
effects.
* Save the version number of `rstan` in element `version` of `brmsfit` objects.


### Bug fixes

* Fix a rare error when predicting `von_mises` models thanks to John Kirwan.



# brms 1.4.0

### Features

* Fit quantile regression models via family `asym_laplace` (asymmetric Laplace
distribution).
* Specify non-linear models in a (hopefully) more intuitive way using
`brmsformula`.
* Fix auxiliary parameters to certain values through `brmsformula`.
* Allow `family` to be specified in `brmsformula`.
* Introduce family `frechet` for modelling strictly positive responses.
* Allow truncation and censoring at the same time.
* Introduce function `prior_` allowing to specify priors using one-sided
formulas or `quote`.
* Pass priors to `Stan` directly without performing any checks by setting `check
= FALSE` in `set_prior`.
* Introduce method `nsamples` to extract the number of posterior samples.
* Export the main formula parsing function `parse_bf`.
* Add more options to customize two-dimensional surface plots created by
`marginal_effects` or `marginal_smooths`.


### Other changes

* Change structure of `brmsformula` objects to be more reliable and easier to
extend.
* Make sure that parameter `nu` never falls below `1` to reduce convergence
problems when using family `student`.
* Deprecate argument `nonlinear`.
* Deprecate family `geometric`.
* Rename `cov_fixed` to `cor_fixed`.
* Make handling of addition terms more transparent by exporting and documenting
related functions.
* Refactor helper functions of the `fitted` method to be easier to extend in the
future.
* Remove many units tests of internal functions and add tests of user-facing
functions instead.
* Import some generics from `nlme` instead of `lme4` to remove dependency on the
latter one.
* Do not apply `structure` to `NULL` anymore to get rid of warnings in R-devel.


### Bug fixes

* Fix problems when fitting smoothing terms with factors as `by` variables
thanks to Milani Chaloupka.
* Fix a bug that could cause some monotonic effects to be ignored in the `Stan`
code thanks to the GitHub user bschneider.
* Make sure that the data of models with only a single observation are
compatible with the generated `Stan` code.
* Handle argument `algorithm` correctly in `update.brmsfit`.
* Fix a bug sometimes causing an error in `marginal_effects` when using family
`wiener` thanks to Andrew Ellis.
* Fix problems in `fitted` when applied to `zero_inflated_beta` models thanks to
Milani Chaloupka.
* Fix minor problems related to the prediction of autocorrelated models.
* Fix a few minor bugs related to the backwards compatibility of multivariate
and related models fitted with `brms` < 1.0.0.




# brms 1.3.1

### Features

* Introduce the auxiliary parameter `disc` ('discrimination') to be used in
ordinal models. By default it is not estimated but fixed to one.
* Create `marginal_effects` plots of two-way interactions of variables that were
not explicitely modeled as interacting.


### Other changes

* Move `rstan` to 'Imports' and `Rcpp` to 'Depends' in order to avoid loading
`rstan` into the global environment automatically.


### Bug fixes

* Fix a bug leading to unexpected errors in some S3 methods when
applied to ordinal models.




# brms 1.3.0

### Features

* Fit error-in-variables models using function `me` in the model formulae.
* Fit multi-membership models using function `mm` in grouping terms.
* Add families `exgaussian` (exponentially modified Gaussian distribution) and
`wiener` (Wiener diffusion model distribution) specifically suited to handle for
response times.
* Add the `lasso` prior as an alternative to the `horseshoe` prior for sparse
models.
* Add the methods `log_posterior`, `nuts_params`, `rhat`, and `neff_ratio` for
`brmsfit` objects to conveniently access quantities used to diagnose sampling
behavior.
* Combine chains in method `as.mcmc` using argument `combine_chains`.
* Estimate the auxiliary parameter `sigma` in models with known standard errors
of the response by setting argument `sigma` to `TRUE` in addition function `se`.
* Allow visualizing two-dimensional smooths with the `marginal_smooths` method.


### Other changes

* Require argument `data` to be explicitely specified in all user facing
functions.
* Refactor the `stanplot` method to use `bayesplot` on the backend.
* Use the `bayesplot` theme as the default in all plotting functions.
* Add the abbreviations `mo` and `cs` to specify monotonic and category specific
effects respectively.
* Rename generated variables in the data.frames returned by `marginal_effects`
to avoid potential naming conflicts.
* Deprecate argument `cluster` and use the native `cores` argument of `rstan`
instead.
* Remove argument `cluster_type` as it is no longer required to apply forking.
* Remove the deprecated `partial` argument.




# brms 1.2.0

### Features

* Add the new family `hurdle_lognormal` specifically suited for zero-inflated
continuous responses.
* Introduce the `pp_check` method to perform various posterior predictive checks
using the `bayesplot` package.
* Introduce the `marginal_smooths` method to better visualize smooth terms.
* Allow varying the scale of global shrinkage parameter of the `horseshoe`
prior.
* Add functions `prior` and `prior_string` as aliases of `set_prior`, the former
allowing to pass arguments without quotes `""` using non-standard evaluation.
* Introduce four new vignettes explaining how to fit non-linear models,
distributional models, phylogenetic models, and monotonic effects respectively.
* Extend the `coef` method to better handle category specific group-level
effects.
* Introduce the `prior_summary` method for `brmsfit` objects to obtain a summary
of prior distributions applied.
* Sample from the prior of the original population-level intercept when
`sample_prior = TRUE` even in models with an internal temporary intercept used
to improve sampling efficiency.
* Introduce methods `posterior_predict`, `predictive_error` and `log_lik` as
(partial) aliases of `predict`, `residuals`, and `logLik` respectively.


### Other changes

* Improve computation of Bayes factors in the `hypothesis` method to be less
influenced by MCMC error.
* Improve documentation of default priors.
* Refactor internal structure of some formula and prior evaluating functions.
This should not have any user visible effects.
* Use the `bayesplot` package as the new backend of `plot.brmsfit`.


### Bug fixes

* Better mimic `mgcv` when parsing smooth terms to make sure all arguments are
correctly handled.
* Avoid an error occurring during the prediction of new data when grouping
factors with only a single factor level were supplied thanks to Tom Wallis.
* Fix `marginal_effects` to consistently produce plots for all covariates in
non-linear models thanks to David Auty.
* Improve the `update` method to better recognize situations where recompliation
of the `Stan` code is necessary thanks to Raphael P.H.
* Allow to correctly `update` the `sample_prior` argument to value `"only"`.
* Fix an unexpected error occurring in many S3 methods when the thinning rate is
not a divisor of the total number of posterior samples thanks to Paul Zerr.



# brms 1.1.0

### Features

* Estimate monotonic group-level effects.
* Estimate category specific group-level effects.
* Allow `t2` smooth terms based on multiple covariates.
* Estimate interval censored data via the addition argument `cens` in the model
formula.
* Allow to compute `residuals` also based on predicted values instead of fitted
values.


### Other changes

* Use the prefix `bcs` in parameter names of category specific effects and the
prefix `bm` in parameter names of monotonic effects (instead of the prefix `b`)
to simplify their identification.
* Ensure full compatibility with `ggplot2` version 2.2.


### Bug fixes

* Fix a bug that could result in incorrect threshold estimates for `cumulative`
and `sratio` models thanks to Peter Congdon.
* Fix a bug that sometimes kept distributional `gamma` models from being
compiled thanks to Tim Beechey.
* Fix a bug causing an error in `predict` and related methods when two-level
factors or logical variables were used as covariates in non-linear models thanks
to Martin Schmettow.
* Fix a bug causing an error when passing lists to additional arguments of
smoothing functions thanks to Wayne Folta.
* Fix a bug causing an error in the `prior_samples` method for models with
multiple group-level terms that refer to the same grouping factor thanks to
Marco Tullio Liuzza.
* Fix a bug sometimes causing an error when calling `marginal_effects` for
weighted models.



# brms 1.0.1
  \subsection{MINOR CHANGES

* Center design matrices inside the Stan code instead of inside `make_standata`.
* Get rid of several warning messages occurring on CRAN.




# brms 1.0.0

This is one of the largest updates of `brms` since its initial release. In
addition to many new features, the multivariate `'trait'` syntax has been
removed from the package as it was confusing for users, required much special
case coding, and was hard to maintain. See `help(brmsformula)` for details of
the formula syntax applied in `brms`.

### Features

* Allow estimating correlations between group-level effects defined across
multiple formulae (e.g., in non-linear models) by specifying IDs in each
grouping term via an extended `lme4` syntax.
* Implement distributional regression models allowing to fully predict auxiliary
parameters of the response distribution. Among many other possibilities, this
can be used to model heterogeneity of variances.
* Zero-inflated and hurdle models do not use multivariate syntax anymore but
instead have special auxiliary parameters named `zi` and `hu` defining
zero-inflation / hurdle probabilities.
* Implement the `von_mises` family to model circular responses.
* Introduce the `brmsfamily` function for convenient specification of `family`
objects.
* Allow predictions of `t2` smoothing terms for new data.
* Feature vectors as arguments for the addition argument `trunc` in order to
model varying truncation points.


### Other changes

* Remove the `cauchy` family after several months of deprecation.
* Make sure that group-level parameter names are unambiguous by adding double
underscores thanks to the idea of the GitHub user schmettow.
* The `predict` method now returns predicted probabilities instead of absolute
frequencies of samples for ordinal and categorical models.
* Compute the linear predictor in the model block of the Stan program instead of
in the transformed parameters block. This avoids saving samples of unnecessary
parameters to disk. Thanks goes to Rick Arrano for pointing me to this issue.
* Colour points in `marginal_effects` plots if sensible.
* Set the default of the `robust` argument to `TRUE` in
`marginal_effects.brmsfit`.


### Bug fixes

* Fix a bug that could occur when predicting factorial response variables for
new data. Only affects categorical and ordinal models.
* Fix a bug that could lead to duplicated variable names in the Stan code when
sampling from priors in non-linear models thanks to Tom Wallis.
* Fix problems when trying to pointwise evaluate non-linear formulae in
`logLik.brmsfit` thanks to Tom Wallis.
* Ensure full compatibility of the `ranef` and `coef` methods with non-linear
models.
* Fix problems that occasionally occurred when handling `dplyr` datasets thanks
to the GitHub user Atan1988.


# brms 0.10.0
### Features

* Add support for generalized additive mixed models (GAMMs). Smoothing terms can
be specified using the `s` and `t2` functions in the model formula.
* Introduce `as.data.frame` and `as.matrix` methods for `brmsfit` objects.

### Other changes

* The `gaussian("log")` family no longer implies a log-normal distribution, but
a normal distribution with log-link to match the behavior of `glm`. The
log-normal distribution can now be specified via family `lognormal`.
* Update syntax of `Stan` models to match the recommended syntax of `Stan` 2.10.


### Bug fixes

* The `ngrps` method should now always return the correct result for non-linear
models.
* Fix problems in `marginal_effects` for models using the reserved variable
`intercept` thanks to Frederik Aust.
* Fix a bug in the `print` method of `brmshypothesis` objects that could lead to
duplicated and thus invalid row names.
* Residual standard deviation parameters of multivariate models are again
correctly displayed in the output of the `summary` method.
* Fix problems when using variational Bayes algorithms with `brms` while having
`rstan` >= 2.10.0 installed thanks to the GitHub user cwerner87.


# brms 0.9.1
### Features

* Allow the '/' symbol in group-level terms in the `formula` argument to
indicate nested grouping structures.
* Allow to compute `WAIC` and `LOO` based on the pointwise log-likelihood using
argument `pointwise` to substantially reduce memory requirements.


### Other changes

* Add horizontal lines to the errorbars in `marginal_effects` plots for factors.


### Bug fixes

* Fix a bug that could lead to a cryptic error message when changing some parts
of the model `formula` using the `update` method.
* Fix a bug that could lead to an error when calling `marginal_effects` for
predictors that were generated with the `base::scale` function thanks to Tom
Wallis.
* Allow interactions of numeric and categorical predictors in `marginal_effects`
to be passed to the `effects` argument in any order.
* Fix a bug that could lead to incorrect results of `predict` and related
methods when called with `newdata` in models using the `poly` function thanks to
Brock Ferguson.
* Make sure that user-specified factor contrasts are always applied in
multivariate models.


# brms 0.9.0
### Features

* Add support for `monotonic` effects allowing to use ordinal predictors without
assuming their categories to be equidistant.
* Apply multivariate formula syntax in categorical models to considerably
increase modeling flexibility.
* Add the addition argument `disp` to define multiplicative factors on
dispersion parameters. For linear models, `disp` applies to the residual
standard deviation `sigma` so that it can be used to weight observations.
* Treat the fixed effects design matrix as sparse by using the `sparse` argument
of `brm`. This can considerably reduce working memory requirements if the
predictors contain many zeros.
* Add the `cor_fixed` correlation structure to allow for fixed user-defined
covariance matrices of the response variable.
* Allow to pass self-defined `Stan` functions via argument `stan_funs` of `brm`.
* Add the `expose_functions` method allowing to expose self-defined `Stan`
functions in `R`.
* Extend the functionality of the `update` method to allow all model parts to be
updated.
* Center the fixed effects design matrix also in multivariate models. This may
lead to increased sampling speed in models with many predictors.


### Other changes

* Refactor `Stan` code and data generating functions to be more consistent and
easier to extent.
* Improve checks of user-define prior specifications.
* Warn about models that have not converged.
* Make sure that regression curves computed by the `marginal_effects` method are
always smooth.
* Allow to define category specific effects in ordinal models directly within
the `formula` argument.


### Bug fixes

* Fix problems in the generated `Stan` code when using very long non-linear
model formulas thanks to Emmanuel Charpentier.
* Fix a bug that prohibited to change priors on single standard deviation
parameters in non-linear models thanks to Emmanuel Charpentier.
* Fix a bug that prohibited to use nested grouping factors in non-linear models
thanks to Tom Wallis.
* Fix a bug in the linear predictor computation within `R`, occurring for ordinal
models with multiple category specific effects. This could lead to incorrect
outputs of `predict`, `fitted`, and `logLik` for these models.
* Make sure that the global `"contrasts"` option is not used when
post-processing a model.


# brms 0.8.0
### Features

* Implement generalized non-linear models, which can be specified with the help
of the `nonlinear` argument in `brm`.
* Compute and plot marginal effects using the `marginal_effects` method thanks
to the help of Ruben Arslan.
* Implement zero-inflated beta models through family `zero_inflated_beta` thanks
to the idea of Ali Roshan Ghias.
* Allow to restrict domain of fixed effects and autocorrelation parameters using
new arguments `lb` and `ub` in function `set_prior` thanks to the idea of Joel
Gombin.
* Add an `as.mcmc` method for compatibility with the `coda` package.
* Allow to call the `WAIC`, `LOO`, and `logLik` methods with new data.


### Other changes

* Make sure that `brms` is fully compatible with `loo` version 0.1.5.
* Optionally define the intercept as an ordinary fixed effect to avoid the
reparametrization via centering of the fixed effects design matrix.
* Do not compute the WAIC in `summary` by default anymore to reduce computation
time of the method for larger models.
* The `cauchy` family is now deprecated and will be removed soon as it often has
convergence issues and not much practical application anyway.
* Change the default settings of the number of chains and warmup samples to the
defaults of `rstan` (i.e., `chains = 4` and `warmup = iter / 2`).
* Do not remove bad behaving chains anymore as they may point to general
convergence problems that are dangerous to ignore.
* Improve flexibility of the `theme` argument in all plotting functions.
* Only show the legend once per page, when computing trace and density plots
with the `plot` method.
* Move code of self-defined `Stan` functions to `inst/chunks` and incorporate
them into the models using `rstan::stanc_builder`. Also, add unit tests for
these functions.


### Bug fixes

* Fix problems when predicting with `newdata` for zero-inflated and hurdle
models thanks to Ruben Arslan.
* Fix problems when predicting with `newdata` if it is a subset of the data
stored in a `brmsfit` object thanks to Ruben Arslan.
* Fix data preparation for multivariate models if some responses are `NA` thanks
to Raphael Royaute.
* Fix a bug in the `predict` method occurring for some multivariate models so
that it now always returns the predictions of all response variables, not just
the first one.
* Fix a bug in the log-likelihood computation of `hurdle_poisson` and
`hurdle_negbinomial` models. This may lead to minor changes in the values
obtained by `WAIC` and `LOO` for these models.
* Fix some backwards compatibility issues of models fitted with version <= 0.5.0
thanks to Ulf Koether.



# brms 0.7.0

### Features

* Use variational inference algorithms as alternative to the NUTS sampler by
specifying argument `algorithm` in the `brm` function.
* Implement beta regression models through family `Beta`.
* Implement zero-inflated binomial models through family
`zero_inflated_binomial`.
* Implement multiplicative effects for family `bernoulli` to fit (among others)
2PL IRT models.
* Generalize the `formula` argument for zero-inflated and hurdle models so that
predictors can be included in only one of the two model parts thanks to the idea
of Wade Blanchard.
* Combine fixed and random effects estimates using the new `coef` method.
* Call the `residuals` method with `newdata` thanks to the idea of Friederike
Holz-Ebeling.
* Allow new levels of random effects grouping factors in the `predict`,
`fitted`, and `residuals` methods using argument `allow_new_levels`.
* Selectively exclude random effects in the `predict`, `fitted`, and `residuals`
methods using argument `re_formula`.
* Add a `plot` method for objects returned by method `hypothesis` to visualize
prior and posterior distributions of the hypotheses being tested.


### Other changes

* Improve evaluation of the response part of the `formula` argument to reliably
allow terms with more than one variable (e.g., `y/x ~ 1`).
* Improve sampling efficiency of models containing many fixed effects through
centering the fixed effects design matrix thanks to Wayne Folta.
* Improve sampling efficiency of models containing uncorrelated random effects
specified by means of `(random || group)` terms in `formula` thanks to Ali
Roshan Ghias.
* Utilize user-defined functions in the `Stan` code of ordinal models to improve
readability as well as sampling efficiency.
* Make sure that model comparisons using `LOO` or `WAIC` are only performed when
models are based on the same responses.
* Use some generic functions of the `lme4` package to avoid unnecessary function
masking. This leads to a change in the argument order of method `VarCorr`.
* Change the `ggplot` theme in the `plot` method through argument `theme`.
* Remove the `n.` prefix in arguments `n.iter`, `n.warmup`, `n.thin`,
`n.chains`, and `n.cluster` of the `brm` function. The old argument names remain
usable as deprecated aliases.
* Amend names of random effects parameters to simplify matching with their
respective grouping factor levels.


### Bug fixes

* Fix a bug in the `hypothesis` method that could cause valid model parameters
to be falsely reported as invalid.
* Fix a bug in the `prior_samples` method that could cause prior samples of
parameters of the same class to be artificially correlated.
* Fix `Stan` code of linear models with moving-average effects and non-identity
link functions so that they no longer contain code related solely to
autoregressive effects.
* Fix a bug in the evaluation of `formula` that could cause complicated random
effects terms to be falsely treated as fixed effects.
* Fix several bugs when calling the `fitted` and `predict` methods with
`newdata` thanks to Ali Roshan Ghias.


# brms 0.6.0

### Features

* Add support for zero-inflated and hurdle models thanks to the idea of Scott
Baldwin.
* Implement inverse gaussian models through family `inverse.gaussian`.
* Allow to specify truncation boundaries of the response variable thanks to the
idea of Maciej Beresewicz.
* Add support for autoregressive (AR) effects of residuals, which can be modeled
using the `cor_ar` and `cor_arma` functions.
* Stationary autoregressive-moving-average (ARMA) effects of order one can now
also be fitted using special covariance matrices.
* Implement multivariate student-t models.
* Binomial and ordinal families now support the `cauchit` link function.
* Allow family functions to be used in the `family` argument.
* Easy access to various `rstan` plotting functions using the `stanplot` method.
* Implement horseshoe priors to model sparsity in fixed effects coefficients
thanks to the idea of Josh Chang.
* Automatically scale default standard deviation priors so that they remain only
weakly informative independent on the response scale.
* Report model weights computed by the `loo` package when comparing multiple
fitted models.


### Other changes

* Separate the fixed effects Intercept from other fixed effects in the `Stan`
code to slightly improve sampling efficiency.
* Move autoregressive (AR) effects of the response from the `cor_ar` to the
`cor_arr` function as the result of implementing AR effects of residuals.
* Improve checks on argument `newdata` used in the `fitted` and `predict`
method.
* Method `standata` is now the only way to extract data that was passed to
`Stan` from a `brmsfit` object.
* Slightly improve `Stan` code for models containing no random effects.
* Change the default prior of the degrees of freedom of the `student` family to
`gamma(2,0.1)`.
* Improve readability of the output of method `VarCorr`.
* Export the `make_stancode` function to give users direct access to `Stan` code
generated by `brms`.
* Rename the `brmdata` function to `make_standata`. The former remains usable as
a deprecated alias.
* Improve documentation to better explain differences in autoregressive effects
across R packages.

### Bug fixes

* Fix a bug that could cause an unexpected error when the `predict` method was
called with `newdata`.
* Avoid side effects of the `rstan` compilation routines that could occasionally
cause R to crash.
* Make `brms` work correctly with `loo` version 0.1.3 thanks to Mauricio Garnier
Villarreal and Jonah Gabry.
* Fix a bug that could cause WAIC and LOO estimates to be slightly incorrect for
`gaussian` models with `log` link.




# brms 0.5.0
### Features

* Compute the Watanabe-Akaike information criterion (WAIC) and leave-one-out
cross-validation (LOO) using the `loo` package.
* Provide an interface to `shinystan` with S3 method `launch_shiny`.
* New functions `get_prior` and `set_prior` to make prior specifications easier.
* Log-likelihood values and posterior predictive samples can now be calculated
within R after the model has been fitted.
* Make predictions based on new data using S3 method `predict`.
* Allow for customized covariance structures of grouping factors with multiple
random effects.
* New S3 methods `fitted` and `residuals` to compute fitted values and
residuals, respectively.


### Other changes

* Arguments `WAIC` and `predict` are removed from the `brm` function, as they
are no longer necessary.
* New argument `cluster_type` in function `brm` allowing to choose the cluster
type created by the parallel package.
* Remove chains that fail to initialize while sampling in parallel leaving the
other chains untouched.
* Redesign trace and density plots to be faster and more stable.
* S3 method `VarCorr` now always returns covariance matrices regardless of
whether correlations were estimated.


### Bug fixes

* Fix a bug in S3 method `hypothesis` related to the calculation of
Bayes-factors for point hypotheses.
* User-defined covariance matrices that are not strictly positive definite for
numerical reasons should now be handled correctly.
* Fix problems when a factor is used as fixed effect and as random effects
grouping variable at the same time thanks to Ulf Koether.
* Fix minor issues with internal parameter naming.
* Perform additional checking on user defined priors.


# brms 0.4.1

### Features

* Allow for sampling from all specified proper priors in the model.
* Compute Bayes-factors for point hypotheses in S3 method `hypothesis`.

### Bug fixes

* Fix a bug that could cause an error for models with multiple grouping factors
thanks to Jonathan Williams.
* Fix a bug that could cause an error for weighted poisson and exponential
models.


# brms 0.4.0

### Features

* Implement the Watanabe-Akaike Information Criterion (WAIC).
* Implement the `||`-syntax for random effects allowing for the estimation of
random effects standard deviations without the estimation of correlations.
* Allow to combine multiple grouping factors within one random effects argument
using the interaction symbol `:`.
* Generalize S3 method `hypothesis` to be used with all parameter classes not
just fixed effects. In addition, one-sided hypothesis testing is now possible.
* Introduce new family `multigaussian` allowing for multivariate normal
regression.
* Introduce new family `bernoulli` for dichotomous response variables as a more
efficient alternative to families `binomial` or `categorical` in this special
case.


### Other changes

* Slightly change the internal structure of brms to reflect that `rstan` is
finally on CRAN.
* Thoroughly check validity of the response variable before the data is passed
to `Stan`.
* Prohibit variable names containing double underscores `__` to avoid naming
conflicts.
* Allow function calls with several arguments (e.g. `poly(x,3)`) in the formula
argument of function `brm`.
* Always center random effects estimates returned by S3 method `ranef` around
zero.
* Prevent the use of customized covariance matrices for grouping factors with
multiple random effects for now.
* Remove any experimental `JAGS` code from the package.


### Bug fixes

* Fix a bug in S3 method `hypothesis` leading to an error when numbers with
decimal places were used in the formulation of the hypotheses.
* Fix a bug in S3 method `ranef` that caused an error for grouping factors with
only one random effect.
* Fix a bug that could cause the fixed intercept to be wrongly estimated in the
presence of multiple random intercepts thanks to Jarrod Hadfield.



# brms 0.3.0

### Features

* Introduce new methods `parnames` and `posterior_samples` for class 'brmsfit'
to extract parameter names and posterior samples for given parameters,
respectively.
* Introduce new method `hypothesis` for class `brmsfit` allowing to test
non-linear hypotheses concerning fixed effects.
* Introduce new argument `addition` in function brm to get a more flexible
approach in specifying additional information on the response variable (e.g.,
standard errors for meta-analysis). Alternatively, this information can also be
passed to the `formula` argument directly.
* Introduce weighted and censored regressions through argument `addition` of
function brm.
* Introduce new argument `cov.ranef` in the `brm` function allowing for
customized covariance structures of random effects thanks to the idea of Boby
Mathew.
* Introduce new argument `autocor` in function brm allowing for autocorrelation
of the response variable.
* Introduce new functions `cor.ar`, `cor.ma`, and `cor.arma`, to be used with
argument `autocor` for modeling autoregressive, moving-average, and
autoregressive-moving-average models.


### Other changes

* Amend parametrization of random effects to increase efficiency of the sampling
algorithms.
* Improve vectorization of sampling statements.


### Bug fixes

* Fix a bug that could cause an error when fitting poisson models while
`predict = TRUE`.
* Fix a bug that caused an error when sampling only one chain while `silent =
TRUE`.


# brms 0.2.0

### Features

* New S3 class `brmsfit` to be returned by the `brm` function.
* New methods for class `brmsfit`: `summary`, `print`, `plot`, `predict`,
`fixef`, `ranef`, `VarCorr`, `nobs`, `ngrps`, and `formula`.
* Introduce new argument `silent` in the `brm` function, allowing to suppress
most of `Stan`'s intermediate output.
* Introduce new families `negbinomial` (negative binomial) and `geometric` to
allow for more flexibility in modeling count data.


### Other changes

* Amend warning and error messages to make them more informative.
* Correct examples in the documentation.
* Extend the README file.


### Bug fixes

* Fix a bug that caused problems when formulas contained more complicated
function calls.
* Fix a bug that caused an error when posterior predictives were sampled for
family `cumulative`.
* Fix a bug that prohibited to use of improper flat priors for parameters that
have proper priors by default.


# brms 0.1.0

* Initial release version


