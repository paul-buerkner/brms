# Extract Posterior Draws

Extract posterior draws in conventional formats as data.frames,
matrices, or arrays.

## Usage

``` r
# S3 method for class 'brmsfit'
as.data.frame(
  x,
  row.names = NULL,
  optional = TRUE,
  pars = NA,
  variable = NULL,
  draw = NULL,
  subset = NULL,
  ...
)

# S3 method for class 'brmsfit'
as.matrix(x, pars = NA, variable = NULL, draw = NULL, subset = NULL, ...)

# S3 method for class 'brmsfit'
as.array(x, pars = NA, variable = NULL, draw = NULL, subset = NULL, ...)
```

## Arguments

- x:

  A `brmsfit` object or another R object for which the methods are
  defined.

- row.names, optional:

  Unused and only added for consistency with the
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) generic.

- pars:

  Deprecated alias of `variable`. For reasons of backwards
  compatibility, `pars` is interpreted as a vector of regular
  expressions by default unless `fixed = TRUE` is specified.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

- draw:

  The draw indices to be select. Subsetting draw indices will lead to an
  automatic merging of chains.

- subset:

  Deprecated alias of `draw`.

- ...:

  Further arguments to be passed to the corresponding
  [`as_draws_*`](https://paulbuerkner.com/brms/reference/draws-brms.md)
  methods as well as to
  [`subset_draws`](https://mc-stan.org/posterior/reference/subset_draws.html).

## Value

A data.frame, matrix, or array containing the posterior draws.

## See also

[`as_draws`](https://paulbuerkner.com/brms/reference/draws-brms.md),
[`subset_draws`](https://mc-stan.org/posterior/reference/subset_draws.html)
