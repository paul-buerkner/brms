# Zero-Inflated Distributions

Density and distribution functions for zero-inflated distributions.

## Usage

``` r
dzero_inflated_poisson(x, lambda, zi, log = FALSE)

pzero_inflated_poisson(q, lambda, zi, lower.tail = TRUE, log.p = FALSE)

dzero_inflated_negbinomial(x, mu, shape, zi, log = FALSE)

pzero_inflated_negbinomial(q, mu, shape, zi, lower.tail = TRUE, log.p = FALSE)

dzero_inflated_binomial(x, size, prob, zi, log = FALSE)

pzero_inflated_binomial(q, size, prob, zi, lower.tail = TRUE, log.p = FALSE)

dzero_inflated_beta_binomial(x, size, mu, phi, zi, log = FALSE)

pzero_inflated_beta_binomial(
  q,
  size,
  mu,
  phi,
  zi,
  lower.tail = TRUE,
  log.p = FALSE
)

dzero_inflated_beta(x, shape1, shape2, zi, log = FALSE)

pzero_inflated_beta(q, shape1, shape2, zi, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- x:

  Vector of quantiles.

- zi:

  zero-inflation probability

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- q:

  Vector of quantiles.

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X
  \> x) .

- log.p:

  Logical; If `TRUE`, values are returned on the log scale.

- mu, lambda:

  location parameter

- shape, shape1, shape2:

  shape parameter

- size:

  number of trials

- prob:

  probability of success on each trial

- phi:

  precision parameter

## Details

The density of a zero-inflated distribution can be specified as follows.
If \\x = 0\\ set \\f(x) = \theta + (1 - \theta) \* g(0)\\. Else set
\\f(x) = (1 - \theta) \* g(x)\\, where \\g(x)\\ is the density of the
non-zero-inflated part.
