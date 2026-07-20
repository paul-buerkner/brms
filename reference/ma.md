# Set up MA(q) correlation structures

Set up a moving average (MA) term of order q in brms. The function does
not evaluate its arguments – it exists purely to help set up a model
with MA terms.

## Usage

``` r
ma(time = NA, gr = NA, q = 1, cov = FALSE)
```

## Arguments

- time:

  An optional time variable specifying the time ordering of the
  observations. By default, the existing order of the observations in
  the data is used.

- gr:

  An optional grouping variable. If specified, the correlation structure
  is assumed to apply only to observations within the same grouping
  level.

- q:

  A non-negative integer specifying the moving average (MA) order of the
  ARMA structure. Default is `1`.

- cov:

  A flag indicating whether ARMA effects should be estimated by means of
  residual covariance matrices. This is currently only possible for
  stationary ARMA effects of order 1. If the model family does not have
  natural residuals, latent residuals are added automatically. If
  `FALSE` (the default), a regression formulation is used that is
  considerably faster and allows for ARMA effects of order higher than 1
  but is only available for `gaussian` models and some of its
  generalizations.

## Value

An object of class `'arma_term'`, which is a list of arguments to be
interpreted by the formula parsing functions of brms.

## See also

[`autocor-terms`](https://paulbuerkner.com/brms/reference/autocor-terms.md),
[`arma`](https://paulbuerkner.com/brms/reference/arma.md),
[`ar`](https://paulbuerkner.com/brms/reference/ar.md)

## Examples

``` r
# \dontrun{
data("LakeHuron")
LakeHuron <- as.data.frame(LakeHuron)
fit <- brm(x ~ ma(p = 2), data = LakeHuron)
#> Error in ma(p = 2): unused argument (p = 2)
summary(fit)
#> Error: object 'fit' not found
# }
```
