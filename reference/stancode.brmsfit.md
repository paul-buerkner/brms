# Extract Stan code from `brmsfit` objects

Extract Stan code from a fitted brms model.

## Usage

``` r
# S3 method for class 'brmsfit'
stancode(
  object,
  version = TRUE,
  regenerate = NULL,
  threads = NULL,
  backend = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

- version:

  Logical; indicates if the first line containing the brms version
  number should be included. Defaults to `TRUE`.

- regenerate:

  Logical; indicates if the Stan code should be regenerated with the
  current brms version. By default, `regenerate` will be `FALSE` unless
  required to be `TRUE` by other arguments.

- threads:

  Controls whether the Stan code should be threaded. See
  [`threading`](https://paulbuerkner.com/brms/reference/threading.md)
  for details.

- backend:

  Controls the Stan backend. See
  [`brm`](https://paulbuerkner.com/brms/reference/brm.md) for details.

- ...:

  Further arguments passed to
  [`stancode`](https://paulbuerkner.com/brms/reference/stancode.default.md)
  if the Stan code is regenerated.

## Value

Stan code for further processing.
