# Extract data passed to Stan from `brmsfit` objects

Extract all data that was used by Stan to fit a brms model.

## Usage

``` r
# S3 method for class 'brmsfit'
standata(
  object,
  newdata = NULL,
  re_formula = NULL,
  newdata2 = NULL,
  new_objects = NULL,
  incl_autocor = TRUE,
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

- newdata2:

  A named `list` of objects containing new data, which cannot be passed
  via argument `newdata`. Required for some objects used in
  autocorrelation structures, or
  [`stanvars`](https://paulbuerkner.com/brms/reference/stanvar.md).

- new_objects:

  Deprecated alias of `newdata2`.

- incl_autocor:

  A flag indicating if correlation structures originally specified via
  `autocor` should be included in the predictions. Defaults to `TRUE`.

- ...:

  More arguments passed to
  [`standata.default`](https://paulbuerkner.com/brms/reference/standata.default.md).
  and
  [`validate_newdata`](https://paulbuerkner.com/brms/reference/validate_newdata.md).

## Value

A named list containing the data passed to Stan.
