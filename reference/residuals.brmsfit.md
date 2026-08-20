# Posterior Draws of Residuals/Predictive Errors

This method is an alias of
[`predictive_error.brmsfit`](https://paulbuerkner.com/brms/reference/predictive_error.brmsfit.md)
with additional arguments for obtaining summaries of the computed draws.

## Usage

``` r
# S3 method for class 'brmsfit'
residuals(
  object,
  newdata = NULL,
  re_formula = NULL,
  method = "posterior_predict",
  type = c("ordinary", "pearson"),
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
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

- method:

  Method used to obtain predictions. Can be set to `"posterior_predict"`
  (the default), `"posterior_epred"`, or `"posterior_linpred"`. For more
  details, see the respective function documentations.

- type:

  The type of the residuals, either `"ordinary"` or `"pearson"`. More
  information is provided under 'Details'.

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

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`..

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead. Only used if `summary` is `TRUE`.

- probs:

  The percentiles to be computed by the `quantile` function. Only used
  if `summary` is `TRUE`.

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

An `array` of predictive error/residual draws. If `summary = FALSE` the
output resembles those of
[`predictive_error.brmsfit`](https://paulbuerkner.com/brms/reference/predictive_error.brmsfit.md).
If `summary = TRUE` the output is an N x E matrix, where N is the number
of observations and E denotes the summary statistics computed from the
draws.

## Details

Residuals of type `'ordinary'` are of the form \\R = Y - Yrep\\, where
\\Y\\ is the observed and \\Yrep\\ is the predicted response. Residuals
of type `pearson` are of the form \\R = (Y - Yrep) / SD(Yrep)\\, where
\\SD(Yrep)\\ is an estimate of the standard deviation of \\Yrep\\.

## Examples

``` r
# \dontrun{
## fit a model
fit <- brm(rating ~ treat + period + carry + (1|subject),
           data = inhaler, cores = 2)
#> Compiling Stan program...
#> Start sampling
#> 

## extract residuals/predictive errors
res <- residuals(fit)
head(res)
# }
```
