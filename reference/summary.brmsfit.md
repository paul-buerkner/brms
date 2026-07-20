# Create a summary of a fitted model represented by a `brmsfit` object

Create a summary of a fitted model represented by a `brmsfit` object

## Usage

``` r
# S3 method for class 'brmsfit'
summary(
  object,
  priors = FALSE,
  prob = 0.95,
  robust = FALSE,
  mc_se = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

- priors:

  Logical; Indicating if priors should be included in the summary.
  Default is `FALSE`.

- prob:

  A value between 0 and 1 indicating the desired probability to be
  covered by the uncertainty intervals. The default is 0.95.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead.

- mc_se:

  Logical; Indicating if the uncertainty in `Estimate` caused by the
  MCMC sampling should be shown in the summary. Defaults to `FALSE`.

- ...:

  Other potential arguments

## Details

The convergence diagnostics `Rhat`, `Bulk_ESS`, and `Tail_ESS` are
described in detail by Vehtari et al. (2020).

## References

Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner (2020). Rank-normalization, folding, and
localization: An improved R-hat for assessing convergence of MCMC.
*Bayesian Analysis*. 1–28. `doi:10.1214/20-BA1221`
