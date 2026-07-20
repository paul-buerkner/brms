# The Shifted Log Normal Distribution

Density, distribution function, quantile function and random generation
for the shifted log normal distribution with mean `meanlog`, standard
deviation `sdlog`, and shift parameter `shift`.

## Usage

``` r
dshifted_lnorm(x, meanlog = 0, sdlog = 1, shift = 0, log = FALSE)

pshifted_lnorm(
  q,
  meanlog = 0,
  sdlog = 1,
  shift = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

qshifted_lnorm(
  p,
  meanlog = 0,
  sdlog = 1,
  shift = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

rshifted_lnorm(n, meanlog = 0, sdlog = 1, shift = 0)
```

## Arguments

- x, q:

  Vector of quantiles.

- meanlog:

  Vector of means.

- sdlog:

  Vector of standard deviations.

- shift:

  Vector of shifts.

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

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.
