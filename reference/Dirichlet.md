# The Dirichlet Distribution

Density function and random number generation for the dirichlet
distribution with shape parameter vector `alpha`.

## Usage

``` r
ddirichlet(x, alpha, log = FALSE)

rdirichlet(n, alpha)
```

## Arguments

- x:

  Matrix of quantiles. Each row corresponds to one probability vector.

- alpha:

  Matrix of positive shape parameters. Each row corresponds to one
  probability vector.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- n:

  Number of draws to sample from the distribution.

## Details

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.
