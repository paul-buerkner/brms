# Compute the Pointwise Log-Likelihood

Compute the Pointwise Log-Likelihood

## Usage

``` r
# S3 method for class 'brmsfit'
log_lik(
  object,
  newdata = NULL,
  re_formula = NULL,
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  pointwise = FALSE,
  combine = TRUE,
  add_point_estimate = FALSE,
  cores = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model object of class `brmsfit`.

- newdata:

  An optional data.frame for which to evaluate predictions. If `NULL`
  (default), the original data of the model is used. `NA` values within
  factors (excluding grouping variables) are interpreted as if all dummy
  variables of this factor are zero. This allows, for instance, to make
  predictions of the grand mean when using sum coding. `NA` values
  within grouping variables are treated as a new level.

- re_formula:

  formula containing group-level effects to be considered in the
  prediction. If `NULL` (default), include all group-level effects; if
  `NA` or `~0`, include no group-level effects.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- ndraws:

  Positive integer indicating how many posterior draws should be used.
  If `NULL` (the default) all draws are used. Ignored if `draw_ids` is
  not `NULL`.

- draw_ids:

  An integer vector specifying the posterior draws to be used. If `NULL`
  (the default), all draws are used.

- pointwise:

  A flag indicating whether to compute the full log-likelihood matrix at
  once (the default), or just return the likelihood function along with
  all data and draws required to compute the log-likelihood separately
  for each observation. The latter option is rarely useful when calling
  `log_lik` directly, but rather when computing
  [`waic`](https://paulbuerkner.com/brms/reference/waic.brmsfit.md) or
  [`loo`](https://paulbuerkner.com/brms/reference/loo.brmsfit.md).

- combine:

  Only relevant in multivariate models. Indicates if the log-likelihoods
  of the submodels should be combined per observation (i.e. added
  together; the default) or if the log-likelihoods should be returned
  separately.

- add_point_estimate:

  For internal use only. Ensures compatibility with the
  [`loo_subsample`](https://paulbuerkner.com/brms/reference/loo_subsample.brmsfit.md)
  method.

- cores:

  Number of cores (defaults to `1`). On non-Windows systems, this
  argument can be set globally via the `mc.cores` option.

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

Usually, an S x N matrix containing the pointwise log-likelihood draws,
where S is the number of draws and N is the number of observations in
the data. For multivariate models and if `combine` is `FALSE`, an S x N
x R array is returned, where R is the number of response variables. If
`pointwise = TRUE`, the output is a function with a `draws` attribute
containing all relevant data and posterior draws.

## Details

`NA` values within factors in `newdata`, are interpreted as if all dummy
variables of this factor are zero. This allows, for instance, to make
predictions of the grand mean when using sum coding.

In multilevel models, it is possible to allow new levels of grouping
factors to be used in the predictions. This can be controlled via
argument `allow_new_levels`. New levels can be sampled in multiple ways,
which can be controlled via argument `sample_new_levels`. Both of these
arguments are documented in
[`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
along with several other useful arguments to control specific aspects of
the predictions.
