# Ordered Beta Distribution

Density, distribution function, quantile function and random generation
for the ordered beta distribution as described in Kubinec (2023).

## Usage

``` r
dordbeta(x, mu, phi, xi, kappa, log = FALSE)

pordbeta(q, mu, phi, xi, kappa, lower.tail = TRUE, log.p = FALSE)

qordbeta(p, mu, phi, xi, kappa, lower.tail = TRUE, log.p = FALSE)

rordbeta(n, mu, phi, xi, kappa)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  mean parameter (on response scale, i.e., in (0, 1))

- phi:

  precision parameter of the beta distribution

- xi:

  cutpoint of boundary at 0 (latent scale)

- kappa:

  offset to xi for cutpoint of boundary at 1 (i.e., cutpoint = xi +
  kappa)

- log, log.p:

  logical; if TRUE, probabilities are given as log(p)

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(X \le x)\\,
  otherwise, \\P(X \> x)\\

- p:

  vector of probabilities

- n:

  number of observations

## Details

The ordered beta distribution combines a beta distribution for values in
(0, 1) with point masses at 0 and 1. The probability of each component
is determined by cutpoint thresholds: P(Y=0) = logit^-1(xi - logit(mu)),
P(Y=1) = 1 - logit^-1((xi + kappa) - logit(mu)), and P(0 \< Y \< 1) is
the remaining probability. For the continuous component, the mean is mu
and precision is phi.

## References

Kubinec, R. (2023). Ordered Beta Regression: A Parsimonious,
Well-Fitting Model for Continuous Data with Lower and Upper Bounds.
*Political Analysis*, 31(4), 519-536.
[doi:10.1017/pan.2022.20](https://doi.org/10.1017/pan.2022.20)
