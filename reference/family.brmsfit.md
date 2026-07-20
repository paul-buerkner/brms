# Extract Model Family Objects

Extract Model Family Objects

## Usage

``` r
# S3 method for class 'brmsfit'
family(object, resp = NULL, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- ...:

  Currently unused.

## Value

A `brmsfamily` object or a list of such objects for multivariate models.
