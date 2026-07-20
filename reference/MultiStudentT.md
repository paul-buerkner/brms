# The Multivariate Student-t Distribution

Density function and random generation for the multivariate Student-t
distribution with location vector `mu`, covariance matrix `Sigma`, and
degrees of freedom `df`.

## Usage

``` r
dmulti_student_t(x, df, mu, Sigma, log = FALSE, check = FALSE)

rmulti_student_t(n, df, mu, Sigma, check = FALSE)
```

## Arguments

- x:

  Vector or matrix of quantiles. If `x` is a matrix, each row is taken
  to be a quantile.

- df:

  Vector of degrees of freedom.

- mu:

  Location vector with length equal to the number of dimensions.

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
