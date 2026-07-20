# Add model fit criteria to model objects

Deprecated aliases of
[`add_criterion`](https://paulbuerkner.com/brms/reference/add_criterion.md).

## Usage

``` r
add_loo(x, model_name = NULL, ...)

add_waic(x, model_name = NULL, ...)

add_ic(x, ...)

# S3 method for class 'brmsfit'
add_ic(x, ic = "loo", model_name = NULL, ...)

add_ic(x, ...) <- value
```

## Arguments

- x:

  An R object typically of class `brmsfit`.

- model_name:

  Optional name of the model. If `NULL` (the default) the name is taken
  from the call to `x`.

- ...:

  Further arguments passed to the underlying functions computing the
  model fit criteria. If you are recomputing an already stored criterion
  with other `...` arguments, make sure to set `overwrite = TRUE`.

- ic, value:

  Names of model fit criteria to compute. Currently supported are
  `"loo"`, `"waic"`, `"kfold"`, `"R2"` (R-squared), and `"marglik"` (log
  marginal likelihood).

## Value

An object of the same class as `x`, but with model fit criteria added
for later usage. Previously computed criterion objects will be
overwritten.
