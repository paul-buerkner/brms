# (Deprecated) AR(p) correlation structure

This function is deprecated. Please see
[`ar`](https://paulbuerkner.com/brms/reference/ar.md) for the new
syntax. This function is a constructor for the `cor_arma` class,
allowing for autoregression terms only.

## Usage

``` r
cor_ar(formula = ~1, p = 1, cov = FALSE)
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
  ARMA structure. Default is 1.

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

An object of class `cor_arma` containing solely autoregression terms.

## Details

AR refers to autoregressive effects of residuals, which is what is
typically understood as autoregressive effects. However, one may also
model autoregressive effects of the response variable, which is called
ARR in brms.

## See also

[`cor_arma`](https://paulbuerkner.com/brms/reference/cor_arma.md)

## Examples

``` r
cor_ar(~visit|patient, p = 2)
#> arma(~visit | patient, 2, 0)
```
