# (Deprecated) Correlation structure classes for the brms package

Classes of correlation structures available in the brms package.
`cor_brms` is not a correlation structure itself, but the class common
to all correlation structures implemented in brms.

## Available correlation structures

- cor_arma:

  autoregressive-moving average (ARMA) structure, with arbitrary orders
  for the autoregressive and moving average components

- cor_ar:

  autoregressive (AR) structure of arbitrary order

- cor_ma:

  moving average (MA) structure of arbitrary order

- cor_car:

  Spatial conditional autoregressive (CAR) structure

- cor_sar:

  Spatial simultaneous autoregressive (SAR) structure

- cor_fixed:

  fixed user-defined covariance structure

## See also

[`cor_arma`](https://paulbuerkner.com/brms/reference/cor_arma.md)`, `[`cor_ar`](https://paulbuerkner.com/brms/reference/cor_ar.md)`, `[`cor_ma`](https://paulbuerkner.com/brms/reference/cor_ma.md)`, `[`cor_car`](https://paulbuerkner.com/brms/reference/cor_car.md)`, `[`cor_sar`](https://paulbuerkner.com/brms/reference/cor_sar.md)`, `[`cor_fixed`](https://paulbuerkner.com/brms/reference/cor_fixed.md)
