# The von Mises Distribution

Density, distribution function, and random generation for the von Mises
distribution with location `mu`, and precision `kappa`.

## Usage

``` r
dvon_mises(x, mu, kappa, log = FALSE)

pvon_mises(q, mu, kappa, lower.tail = TRUE, log.p = FALSE, acc = 1e-20)

rvon_mises(n, mu, kappa)
```

## Arguments

- x, q:

  Vector of quantiles between `-pi` and `pi`.

- mu:

  Vector of location values.

- kappa:

  Vector of precision values.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X
  \> x) .

- log.p:

  Logical; If `TRUE`, values are returned on the log scale.

- acc:

  Accuracy of numerical approximations.

- n:

  Number of draws to sample from the distribution.

## Details

See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.
