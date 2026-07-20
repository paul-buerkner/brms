# The Exponentially Modified Gaussian Distribution

Density, distribution function, and random generation for the
exponentially modified Gaussian distribution with mean `mu` and standard
deviation `sigma` of the gaussian component, as well as scale `beta` of
the exponential component.

## Usage

``` r
dexgaussian(x, mu, sigma, beta, log = FALSE)

pexgaussian(q, mu, sigma, beta, lower.tail = TRUE, log.p = FALSE)

rexgaussian(n, mu, sigma, beta)
```

## Arguments

- x, q:

  Vector of quantiles.

- mu:

  Vector of means of the combined distribution.

- sigma:

  Vector of standard deviations of the gaussian component.

- beta:

  Vector of scales of the exponential component.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X
  \> x) .

- log.p:

  Logical; If `TRUE`, values are returned on the log scale.

- n:

  Number of draws to sample from the distribution.

## Details

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.
