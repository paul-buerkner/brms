# Validate New Data

Validate new data passed to post-processing methods of brms. Unless you
are a package developer, you will rarely need to call `validate_newdata`
directly.

## Usage

``` r
validate_newdata(
  newdata,
  object,
  re_formula = NULL,
  allow_new_levels = FALSE,
  newdata2 = NULL,
  resp = NULL,
  check_response = TRUE,
  incl_autocor = TRUE,
  group_vars = NULL,
  req_vars = NULL,
  ...
)
```

## Arguments

- newdata:

  A `data.frame` containing new data to be validated.

- object:

  A `brmsfit` object.

- re_formula:

  formula containing group-level effects to be considered in the
  prediction. If `NULL` (default), include all group-level effects; if
  `NA` or `~0`, include no group-level effects.

- allow_new_levels:

  A flag indicating if new levels of group-level effects are allowed
  (defaults to `FALSE`). Only relevant if `newdata` is provided.

- newdata2:

  A named `list` of objects containing new data, which cannot be passed
  via argument `newdata`. Required for some objects used in
  autocorrelation structures, or
  [`stanvars`](https://paulbuerkner.com/brms/reference/stanvar.md).

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- check_response:

  Logical; Indicates if response variables should be checked as well.
  Defaults to `TRUE`.

- incl_autocor:

  A flag indicating if correlation structures originally specified via
  `autocor` should be included in the predictions. Defaults to `TRUE`.

- group_vars:

  Optional names of grouping variables to be validated. Defaults to all
  grouping variables in the model.

- req_vars:

  Optional names of variables required in `newdata`. If `NULL` (the
  default), all variables in the original data are required (unless
  ignored for some other reason).

- ...:

  Currently ignored.

## Value

A validated `'data.frame'` based on `newdata`.
