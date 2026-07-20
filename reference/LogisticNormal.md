# The (Multivariate) Logistic Normal Distribution

Density function and random generation for the (multivariate) logistic
normal distribution with latent mean vector `mu` and covariance matrix
`Sigma`.

## Usage

``` r
dlogistic_normal(x, mu, Sigma, refcat = 1, log = FALSE, check = FALSE)

rlogistic_normal(n, mu, Sigma, refcat = 1, check = FALSE)
```

## Arguments

- x:

  Vector or matrix of quantiles. If `x` is a matrix, each row is taken
  to be a quantile.

- mu:

  Mean vector with length equal to the number of dimensions.

- Sigma:

  Covariance matrix.

- refcat:

  A single integer indicating the reference category. Defaults to `1`.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- check:

  Logical; Indicates whether several input checks should be performed.
  Defaults to `FALSE` to improve efficiency.

- n:

  Number of draws to sample from the distribution.
