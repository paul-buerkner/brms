# Posterior Draws of Predictive Errors

Compute posterior draws of predictive errors, that is, observed minus
predicted responses. Can be performed for the data used to fit the model
(posterior predictive checks) or for new data.

## Usage

``` r
# S3 method for class 'brmsfit'
predictive_error(
  object,
  newdata = NULL,
  re_formula = NULL,
  re.form = NULL,
  method = "posterior_predict",
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

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

- re.form:

  Alias of `re_formula`.

- method:

  Method used to obtain predictions. Can be set to `"posterior_predict"`
  (the default), `"posterior_epred"`, or `"posterior_linpred"`. For more
  details, see the respective function documentations.

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

- sort:

  Logical. Only relevant for time series models. Indicating whether to
  return predicted values in the original order (`FALSE`; default) or in
  the order of the time series (`TRUE`).

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

An S x N `array` of predictive error draws, where S is the number of
posterior draws and N is the number of observations.

## Examples

``` r
# \dontrun{
## fit a model
fit <- brm(rating ~ treat + period + carry + (1|subject),
           data = inhaler, cores = 2)
#> Compiling Stan program...
#> Start sampling
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

## extract predictive errors
pe <- predictive_error(fit)
str(pe)
# }
```
