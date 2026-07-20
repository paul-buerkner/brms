# Add compiled rstan models to `brmsfit` objects

Compile a
[`stanmodel`](https://mc-stan.org/rstan/reference/stanmodel-class.html)
and add it to a `brmsfit` object. This enables some advanced
functionality of rstan, most notably
[`log_prob`](https://mc-stan.org/rstan/reference/stanfit-method-logprob.html)
and friends, to be used with brms models fitted with other Stan
backends.

## Usage

``` r
add_rstan_model(x, overwrite = FALSE)
```

## Arguments

- x:

  A `brmsfit` object to be updated.

- overwrite:

  Logical. If `TRUE`, overwrite any existing
  [`stanmodel`](https://mc-stan.org/rstan/reference/stanmodel-class.html).
  Defaults to `FALSE`.

## Value

A (possibly updated) `brmsfit` object.
