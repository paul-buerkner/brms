# (Deprecated) Compound Symmetry (COSY) Correlation Structure

This function is deprecated. Please see
[`cosy`](https://paulbuerkner.com/brms/reference/cosy.md) for the new
syntax. This functions is a constructor for the `cor_cosy` class,
representing a compound symmetry structure corresponding to uniform
correlation.

## Usage

``` r
cor_cosy(formula = ~1)
```

## Arguments

- formula:

  A one sided formula of the form `~ t`, or `~ t | g`, specifying a time
  covariate `t` and, optionally, a grouping factor `g`. A covariate for
  this correlation structure must be integer valued. When a grouping
  factor is present in `formula`, the correlation structure is assumed
  to apply only to observations within the same grouping level;
  observations with different grouping levels are assumed to be
  uncorrelated. Defaults to `~ 1`, which corresponds to using the order
  of the observations in the data as a covariate, and no groups.

## Value

An object of class `cor_cosy`, representing a compound symmetry
correlation structure.

## Examples

``` r
cor_cosy(~ visit | patient)
#> cosy(~visit | patient)
```
