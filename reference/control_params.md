# Extract Control Parameters of the NUTS Sampler

Extract control parameters of the NUTS sampler such as `adapt_delta` or
`max_treedepth`.

## Usage

``` r
control_params(x, ...)

# S3 method for class 'brmsfit'
control_params(x, pars = NULL, ...)
```

## Arguments

- x:

  An R object

- ...:

  Currently ignored.

- pars:

  Optional names of the control parameters to be returned. If `NULL`
  (the default) all control parameters are returned. See
  [`stan`](https://mc-stan.org/rstan/reference/stan.html) for more
  details.

## Value

A named `list` with control parameter values.
