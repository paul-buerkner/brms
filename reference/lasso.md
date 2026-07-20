# (Defunct) Set up a lasso prior in brms

This functionality is no longer supported as of brms version 2.19.2.
Please use the
[`horseshoe`](https://paulbuerkner.com/brms/reference/horseshoe.md) or
[`R2D2`](https://paulbuerkner.com/brms/reference/R2D2.md) shrinkage
priors instead.

## Usage

``` r
lasso(df = 1, scale = 1)
```

## Arguments

- df:

  Degrees of freedom of the chi-square prior of the inverse tuning
  parameter. Defaults to `1`.

- scale:

  Scale of the lasso prior. Defaults to `1`.

## Value

An error indicating that the lasso prior is no longer supported.

## References

Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the
American Statistical Association, 103(482), 681-686.

## See also

[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md),
[`horseshoe`](https://paulbuerkner.com/brms/reference/horseshoe.md),
[`R2D2`](https://paulbuerkner.com/brms/reference/R2D2.md)
