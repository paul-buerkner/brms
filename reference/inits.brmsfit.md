# Extract Initial Values Used for Each Chain

Extract the initial values used by Stan for each chain. This is
currently only available if the model was fitted with the `rstan`
backend.

## Usage

``` r
# S3 method for class 'brmsfit'
inits(x, ...)

inits(x, ...)
```

## Arguments

- x:

  A `brmsfit` object.

- ...:

  Currently ignored.

## Value

The initial values (either user-specified or generated randomly) for all
chains. This is a list with one component per chain. Each component is a
named list containing the initial values for each parameter for the
corresponding chain.
