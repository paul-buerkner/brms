# The Generalized Extreme Value Distribution

Density, distribution function, and random generation for the
generalized extreme value distribution with location `mu`, scale `sigma`
and shape `xi`.

## Usage

``` r
dgen_extreme_value(x, mu = 0, sigma = 1, xi = 0, log = FALSE)

pgen_extreme_value(
  q,
  mu = 0,
  sigma = 1,
  xi = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

qgen_extreme_value(
  p,
  mu = 0,
  sigma = 1,
  xi = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

rgen_extreme_value(n, mu = 0, sigma = 1, xi = 0)
```

## Arguments

- x, q:

  Vector of quantiles.

- mu:

  Vector of locations.

- sigma:

  Vector of scales.

- xi:

  Vector of shapes.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X
  \> x) .

- log.p:

  Logical; If `TRUE`, values are returned on the log scale.

- p:

  Vector of probabilities.

- n:

  Number of draws to sample from the distribution.

## Details

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.
