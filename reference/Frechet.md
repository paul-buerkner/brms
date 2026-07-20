# The Frechet Distribution

Density, distribution function, quantile function and random generation
for the Frechet distribution with location `loc`, scale `scale`, and
shape `shape`.

## Usage

``` r
dfrechet(x, loc = 0, scale = 1, shape = 1, log = FALSE)

pfrechet(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

qfrechet(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

rfrechet(n, loc = 0, scale = 1, shape = 1)
```

## Arguments

- x, q:

  Vector of quantiles.

- loc:

  Vector of locations.

- scale:

  Vector of scales.

- shape:

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
