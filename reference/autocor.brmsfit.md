# (Deprecated) Extract Autocorrelation Objects

(Deprecated) Extract Autocorrelation Objects

## Usage

``` r
# S3 method for class 'brmsfit'
autocor(object, resp = NULL, ...)

autocor(object, ...)
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

A `cor_brms` object or a list of such objects for multivariate models.
Not supported for models fitted with brms 2.11.1 or higher.
