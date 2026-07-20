# (Deprecated) Extract posterior samples for use with the coda package

The `as.mcmc` method is deprecated. We recommend using the more modern
and consistent
[`as_draws_*`](https://paulbuerkner.com/brms/reference/draws-brms.md)
extractor functions of the posterior package instead.

## Usage

``` r
# S3 method for class 'brmsfit'
as.mcmc(
  x,
  pars = NA,
  fixed = FALSE,
  combine_chains = FALSE,
  inc_warmup = FALSE,
  ...
)
```

## Arguments

- x:

  An `R` object typically of class `brmsfit`

- pars:

  Names of parameters for which posterior samples should be returned, as
  given by a character vector or regular expressions. By default, all
  posterior samples of all parameters are extracted.

- fixed:

  Indicates whether parameter names should be matched exactly (`TRUE`)
  or treated as regular expressions (`FALSE`). Default is `FALSE`.

- combine_chains:

  Indicates whether chains should be combined.

- inc_warmup:

  Indicates if the warmup samples should be included. Default is
  `FALSE`. Warmup samples are used to tune the parameters of the
  sampling algorithm and should not be analyzed.

- ...:

  currently unused

## Value

If `combine_chains = TRUE` an `mcmc` object is returned. If
`combine_chains = FALSE` an `mcmc.list` object is returned.
