# (Deprecated) Number of Posterior Samples

Extract the number of posterior samples (draws) stored in a fitted
Bayesian model. Method `nsamples` is deprecated. Please use `ndraws`
instead.

## Usage

``` r
# S3 method for class 'brmsfit'
nsamples(object, subset = NULL, incl_warmup = FALSE, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- subset:

  An optional integer vector defining a subset of samples to be
  considered.

- incl_warmup:

  A flag indicating whether to also count warmup / burn-in samples.

- ...:

  Currently ignored.
