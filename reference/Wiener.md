# The Wiener Diffusion Model Distribution

Density function and random generation for the Wiener diffusion model
distribution with boundary separation `alpha`, non-decision time `tau`,
bias `beta` and drift rate `delta`.

## Usage

``` r
dwiener(
  x,
  alpha,
  tau,
  beta,
  delta,
  resp = 1,
  log = FALSE,
  backend = getOption("wiener_backend", "Rwiener")
)

rwiener(
  n,
  alpha,
  tau,
  beta,
  delta,
  types = c("q", "resp"),
  backend = getOption("wiener_backend", "Rwiener")
)
```

## Arguments

- x:

  Vector of quantiles.

- alpha:

  Boundary separation parameter.

- tau:

  Non-decision time parameter.

- beta:

  Bias parameter.

- delta:

  Drift rate parameter.

- resp:

  Response: `"upper"` or `"lower"`. If no character vector, it is
  coerced to logical where `TRUE` indicates `"upper"` and `FALSE`
  indicates `"lower"`.

- log:

  Logical; If `TRUE`, values are returned on the log scale.

- backend:

  Name of the package to use as backend for the computations. Either
  `"Rwiener"` (the default) or `"rtdists"`. Can be set globally for the
  current R session via the `"wiener_backend"` option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- n:

  Number of draws to sample from the distribution.

- types:

  Which types of responses to return? By default, return both the
  response times `"q"` and the dichotomous responses `"resp"`. If either
  `"q"` or `"resp"`, return only one of the two types.

## Details

These are wrappers around functions of the RWiener or rtdists package
(depending on the chosen `backend`). See
[`vignette("brms_families")`](https://paulbuerkner.com/brms/articles/brms_families.md)
for details on the parameterization.

## See also

`wienerdist`,
[`Diffusion`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
