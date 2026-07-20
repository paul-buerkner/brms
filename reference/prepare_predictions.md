# Prepare Predictions

This method helps in preparing brms models for certain post-processing
tasks most notably various forms of predictions. Unless you are a
package developer, you will rarely need to call `prepare_predictions`
directly.

## Usage

``` r
# S3 method for class 'brmsfit'
prepare_predictions(
  x,
  newdata = NULL,
  re_formula = NULL,
  allow_new_levels = FALSE,
  sample_new_levels = "uncertainty",
  incl_autocor = TRUE,
  oos = NULL,
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  nsamples = NULL,
  subset = NULL,
  nug = NULL,
  smooths_only = FALSE,
  offset = TRUE,
  newdata2 = NULL,
  new_objects = NULL,
  point_estimate = NULL,
  ndraws_point_estimate = 1,
  ...
)

prepare_predictions(x, ...)
```

## Arguments

- x:

  An R object typically of class `'brmsfit'`.

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

- allow_new_levels:

  A flag indicating if new levels of group-level effects are allowed
  (defaults to `FALSE`). Only relevant if `newdata` is provided.

- sample_new_levels:

  Indicates how to sample new levels for grouping factors specified in
  `re_formula`. This argument is only relevant if `newdata` is provided
  and `allow_new_levels` is set to `TRUE`. If `"uncertainty"` (default),
  each posterior sample for a new level is drawn from the posterior
  draws of a randomly chosen existing level. Each posterior sample for a
  new level may be drawn from a different existing level such that the
  resulting set of new posterior draws represents the variation across
  existing levels. If `"gaussian"`, sample new levels from the
  (multivariate) normal distribution implied by the group-level standard
  deviations and correlations. This options may be useful for conducting
  Bayesian power analysis or predicting new levels in situations where
  relatively few levels where observed in the old_data. If
  `"old_levels"`, directly sample new levels from the existing levels,
  where a new level is assigned all of the posterior draws of the same
  (randomly chosen) existing level.

- incl_autocor:

  A flag indicating if correlation structures originally specified via
  `autocor` should be included in the predictions. Defaults to `TRUE`.

- oos:

  Optional indices of observations for which to compute out-of-sample
  rather than in-sample predictions. Only required in models that make
  use of response values to make predictions, that is, currently only
  ARMA models.

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

- nsamples:

  Deprecated alias of `ndraws`.

- subset:

  Deprecated alias of `draw_ids`.

- nug:

  Small positive number for Gaussian process terms only. For numerical
  reasons, the covariance matrix of a Gaussian process might not be
  positive definite. Adding a very small number to the matrix's diagonal
  often solves this problem. If `NULL` (the default), `nug` is chosen
  internally.

- smooths_only:

  Logical; If `TRUE` only predictions related to smoothing splines
  (i.e., `s` or `t2`) will be computed. Defaults to `FALSE`.

- offset:

  Logical; Indicates if offsets should be included in the predictions.
  Defaults to `TRUE`.

- newdata2:

  A named `list` of objects containing new data, which cannot be passed
  via argument `newdata`. Required for some objects used in
  autocorrelation structures, or
  [`stanvars`](https://paulbuerkner.com/brms/reference/stanvar.md).

- new_objects:

  Deprecated alias of `newdata2`.

- point_estimate:

  Shall the returned object contain only point estimates of the
  parameters instead of their posterior draws? Defaults to `NULL` in
  which case no point estimate is computed. Alternatively, may be set to
  `"mean"` or `"median"`. This argument is primarily implemented to
  ensure compatibility with the
  [`loo_subsample`](https://paulbuerkner.com/brms/reference/loo_subsample.brmsfit.md)
  method.

- ndraws_point_estimate:

  Only used if `point_estimate` is not `NULL`. How often shall the point
  estimate's value be repeated? Defaults to `1`.

- ...:

  Further arguments passed to
  [`validate_newdata`](https://paulbuerkner.com/brms/reference/validate_newdata.md).

## Value

An object of class `'brmsprep'` or `'mvbrmsprep'`, depending on whether
a univariate or multivariate model is passed.
