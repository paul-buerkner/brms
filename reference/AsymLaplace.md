# The Asymmetric Laplace Distribution

Density, distribution function, quantile function and random generation
for the asymmetric Laplace distribution with location `mu`, scale
`sigma` and asymmetry parameter `quantile`.

## Usage

``` r
dasym_laplace(x, mu = 0, sigma = 1, quantile = 0.5, log = FALSE)

pasym_laplace(
  q,
  mu = 0,
  sigma = 1,
  quantile = 0.5,
  lower.tail = TRUE,
  log.p = FALSE
)

qasym_laplace(
  p,
  mu = 0,
  sigma = 1,
  quantile = 0.5,
  lower.tail = TRUE,
  log.p = FALSE
)

rasym_laplace(n, mu = 0, sigma = 1, quantile = 0.5)
```

## Arguments

- x, q:

  Vector of quantiles.

- mu:

  Vector of locations.

- sigma:

  Vector of scales.

- quantile:

  Asymmetry parameter corresponding to quantiles in quantile regression
  (hence the name).

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
