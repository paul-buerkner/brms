# (Deprecated) ARMA(p,q) correlation structure

This function is deprecated. Please see
[`arma`](https://paulbuerkner.com/brms/reference/arma.md) for the new
syntax. This functions is a constructor for the `cor_arma` class,
representing an autoregression-moving average correlation structure of
order (p, q).

## Usage

``` r
cor_arma(formula = ~1, p = 0, q = 0, r = 0, cov = FALSE)
```

## Arguments

- formula:

  A one sided formula of the form `~ t`, or `~ t | g`, specifying a time
  covariate `t` and, optionally, a grouping factor `g`. A covariate for
  this correlation structure must be integer valued. When a grouping
  factor is present in `formula`, the correlation structure is assumed
  to apply only to observations within the same grouping level;
  observations with different grouping levels are assumed to be
  uncorrelated. Defaults to `~ 1`, which corresponds to using the order
  of the observations in the data as a covariate, and no groups.

- p:

  A non-negative integer specifying the autoregressive (AR) order of the
  ARMA structure. Default is 0.

- q:

  A non-negative integer specifying the moving average (MA) order of the
  ARMA structure. Default is 0.

- r:

  No longer supported.

- cov:

  A flag indicating whether ARMA effects should be estimated by means of
  residual covariance matrices. This is currently only possible for
  stationary ARMA effects of order 1. If the model family does not have
  natural residuals, latent residuals are added automatically. If
  `FALSE` (the default) a regression formulation is used that is
  considerably faster and allows for ARMA effects of order higher than 1
  but is only available for `gaussian` models and some of its
  generalizations.

## Value

An object of class `cor_arma`, representing an
autoregression-moving-average correlation structure.

## See also

[`cor_ar`](https://paulbuerkner.com/brms/reference/cor_ar.md),
[`cor_ma`](https://paulbuerkner.com/brms/reference/cor_ma.md)

## Examples

``` r
cor_arma(~ visit | patient, p = 2, q = 2)
#> arma(~visit | patient, 2, 2)
```
