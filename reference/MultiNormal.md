# The Multivariate Normal Distribution

Density function and random generation for the multivariate normal
distribution with mean vector `mu` and covariance matrix `Sigma`.

## Usage

``` r
dmulti_normal(x, mu, Sigma, log = FALSE, check = FALSE)

rmulti_normal(n, mu, Sigma, check = FALSE)
```

## Arguments

- x:

  Vector or matrix of quantiles. If `x` is a matrix, each row is taken
  to be a quantile.

- mu:

  Mean vector with length equal to the number of dimensions.

- Sigma:

  Covariance matrix.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- check:

  Logical; Indicates whether several input checks should be performed.
  Defaults to `FALSE` to improve efficiency.

- n:

  Number of draws to sample from the distribution.

## Details

See the Stan user's manual <https://mc-stan.org/docs/> for details on
the parameterization
