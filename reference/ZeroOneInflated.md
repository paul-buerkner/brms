# Zero-One-Inflated Beta Distribution

Density, distribution function and quantile function for the
zero-one-inflated beta distribution.

## Usage

``` r
dzero_one_inflated_beta(x, shape1, shape2, zoi, coi, log = FALSE)

pzero_one_inflated_beta(
  q,
  shape1,
  shape2,
  zoi,
  coi,
  lower.tail = TRUE,
  log.p = FALSE
)

qzero_one_inflated_beta(
  p,
  shape1,
  shape2,
  zoi,
  coi,
  lower.tail = TRUE,
  log.p = FALSE
)
```

## Arguments

- x:

  Vector of quantiles.

- shape1, shape2:

  shape parameters of the beta component

- zoi:

  zero-one-inflation probability

- coi:

  conditional one-inflation probability

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- q:

  Vector of quantiles.

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X
  \> x) .

- log.p:

  Logical; If `TRUE`, values are returned on the log scale.

- p:

  Vector of probabilities.

## Details

With probability \\\zeta\\ the response is 0 or 1, and conditionally on
inflation it is 1 with probability \\\gamma\\. With probability \\1 -
\zeta\\ the response follows a beta distribution.
