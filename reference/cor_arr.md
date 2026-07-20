# (Defunct) ARR correlation structure

The ARR correlation structure is no longer supported.

## Usage

``` r
cor_arr(formula = ~1, r = 1)
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

- r:

  No longer supported.
