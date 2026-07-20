# The Student-t Distribution

Density, distribution function, quantile function and random generation
for the Student-t distribution with location `mu`, scale `sigma`, and
degrees of freedom `df`.

## Usage

``` r
dstudent_t(x, df, mu = 0, sigma = 1, log = FALSE)

pstudent_t(q, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)

qstudent_t(p, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)

rstudent_t(n, df, mu = 0, sigma = 1)
```

## Arguments

- x:

  Vector of quantiles.

- df:

  Vector of degrees of freedom.

- mu:

  Vector of location values.

- sigma:

  Vector of scale values.

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

- n:

  Number of draws to sample from the distribution.

## Details

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.

## See also

[`TDist`](https://rdrr.io/r/stats/TDist.html)
