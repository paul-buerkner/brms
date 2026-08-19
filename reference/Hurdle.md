# Hurdle Distributions

Density, distribution function, quantile function and random generation
for hurdle distributions.

## Usage

``` r
dhurdle_poisson(x, lambda, hu, log = FALSE)

phurdle_poisson(q, lambda, hu, lower.tail = TRUE, log.p = FALSE)

qhurdle_poisson(p, lambda, hu, lower.tail = TRUE, log.p = FALSE)

rhurdle_poisson(n, lambda, hu)

dhurdle_negbinomial(x, mu, shape, hu, log = FALSE)

phurdle_negbinomial(q, mu, shape, hu, lower.tail = TRUE, log.p = FALSE)

qhurdle_negbinomial(p, mu, shape, hu, lower.tail = TRUE, log.p = FALSE)

rhurdle_negbinomial(n, mu, shape, hu)

dhurdle_gamma(x, shape, scale, hu, log = FALSE)

phurdle_gamma(q, shape, scale, hu, lower.tail = TRUE, log.p = FALSE)

qhurdle_gamma(p, shape, scale, hu, lower.tail = TRUE, log.p = FALSE)

rhurdle_gamma(n, shape, scale, hu)

dhurdle_lognormal(x, mu, sigma, hu, log = FALSE)

phurdle_lognormal(q, mu, sigma, hu, lower.tail = TRUE, log.p = FALSE)

qhurdle_lognormal(p, mu, sigma, hu, lower.tail = TRUE, log.p = FALSE)

rhurdle_lognormal(n, mu, sigma, hu)
```

## Arguments

- x:

  Vector of quantiles.

- hu:

  hurdle probability

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

- mu, lambda:

  location parameter

- shape:

  shape parameter

- sigma, scale:

  scale parameter

## Details

The density of a hurdle distribution can be specified as follows. If \\x
= 0\\ set \\f(x) = \theta\\. Else set \\f(x) = (1 - \theta) \* g(x) /
(1 - G(0))\\ where \\g(x)\\ and \\G(x)\\ are the density and
distribution function of the non-hurdle part, respectively.
