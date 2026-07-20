# Restructure Old `brmsfit` Objects

Restructure old `brmsfit` objects to work with the latest brms version.
This function is called internally when applying post-processing
methods. However, in order to avoid unnecessary run time caused by the
restructuring, I recommend explicitly calling `restructure` once per
model after updating brms.

## Usage

``` r
# S3 method for class 'brmsfit'
restructure(x, ...)
```

## Arguments

- x:

  An object of class `brmsfit`.

- ...:

  Currently ignored.

## Value

A `brmsfit` object compatible with the latest version of brms.

## Details

If you are restructuring an old spline model (fitted with brms \<
2.19.3) to avoid prediction inconsistencies between machines (see GitHub
issue \#1465), please make sure to `restructure` your model on the
machine on which it was originally fitted.
