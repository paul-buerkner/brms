# The Beta-binomial Distribution

Cumulative density & mass functions, and random number generation for
the Beta-binomial distribution using the following re-parameterisation
of the [Stan Beta-binomial
definition](https://mc-stan.org/docs/2_29/functions-reference/beta-binomial-distribution.html):

- `mu = alpha * beta` mean probability of trial success.

- `phi = (1 - mu) * beta` precision or over-dispersion, component.

## Usage

``` r
dbeta_binomial(x, size, mu, phi, log = FALSE)

pbeta_binomial(q, size, mu, phi, lower.tail = TRUE, log.p = FALSE)

qbeta_binomial(p, size, mu, phi, lower.tail = TRUE, log.p = FALSE)

rbeta_binomial(n, size, mu, phi)
```

## Arguments

- x, q:

  Vector of quantiles.

- size:

  Vector of number of trials (zero or more).

- mu:

  Vector of means.

- phi:

  Vector of precisions.

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

The quantile function has no known closed form for this parameterization
and is therefore computed numerically via inversion of the cumulative
distribution function over the finite support `0:size`.
