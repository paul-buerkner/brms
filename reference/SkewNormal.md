# The Skew-Normal Distribution

Density, distribution function, and random generation for the
skew-normal distribution with mean `mu`, standard deviation `sigma`, and
skewness `alpha`.

## Usage

``` r
dskew_normal(
  x,
  mu = 0,
  sigma = 1,
  alpha = 0,
  xi = NULL,
  omega = NULL,
  log = FALSE
)

pskew_normal(
  q,
  mu = 0,
  sigma = 1,
  alpha = 0,
  xi = NULL,
  omega = NULL,
  lower.tail = TRUE,
  log.p = FALSE
)

qskew_normal(
  p,
  mu = 0,
  sigma = 1,
  alpha = 0,
  xi = NULL,
  omega = NULL,
  lower.tail = TRUE,
  log.p = FALSE,
  tol = 1e-08
)

rskew_normal(n, mu = 0, sigma = 1, alpha = 0, xi = NULL, omega = NULL)
```

## Arguments

- x, q:

  Vector of quantiles.

- mu:

  Vector of mean values.

- sigma:

  Vector of standard deviation values.

- alpha:

  Vector of skewness values.

- xi:

  Optional vector of location values. If `NULL` (the default), will be
  computed internally.

- omega:

  Optional vector of scale values. If `NULL` (the default), will be
  computed internally.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X
  \> x) .

- log.p:

  Logical; If `TRUE`, values are returned on the log scale.

- p:

  Vector of probabilities.

- tol:

  Tolerance of the approximation used in the computation of quantiles.

- n:

  Number of draws to sample from the distribution.

## Details

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.
