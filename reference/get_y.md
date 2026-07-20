# Extract response values

Extract response values from a
[`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.md)
object.

## Usage

``` r
get_y(x, resp = NULL, sort = FALSE, warn = FALSE, ...)
```

## Arguments

- x:

  A
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.md)
  object.

- resp:

  Optional names of response variables for which to extract values.

- sort:

  Logical. Only relevant for time series models. Indicating whether to
  return predicted values in the original order (`FALSE`; default) or in
  the order of the time series (`TRUE`).

- warn:

  For internal use only.

- ...:

  Further arguments passed to
  [`standata`](https://paulbuerkner.com/brms/reference/standata.md).

## Value

Returns a vector of response values for univariate models and a matrix
of response values with one column per response variable for
multivariate models.
