# Expose user-defined Stan functions

Export user-defined Stan function and optionally vectorize them. For
more details see
[`expose_stan_functions`](https://mc-stan.org/rstan/reference/expose_stan_functions.html).

## Usage

``` r
# S3 method for class 'brmsfit'
expose_functions(x, vectorize = FALSE, env = globalenv(), ...)

expose_functions(x, ...)
```

## Arguments

- x:

  An object of class `brmsfit`.

- vectorize:

  Logical; Indicates if the exposed functions should be vectorized via
  [`Vectorize`](https://rdrr.io/r/base/Vectorize.html). Defaults to
  `FALSE`.

- env:

  Environment where the functions should be made available. Defaults to
  the global environment.

- ...:

  Further arguments passed to
  [`expose_stan_functions`](https://mc-stan.org/rstan/reference/expose_stan_functions.html).
