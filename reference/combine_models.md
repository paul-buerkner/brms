# Combine Models fitted with brms

Combine multiple `brmsfit` objects, which fitted the same model. This is
usefully for instance when having manually run models in parallel.

## Usage

``` r
combine_models(..., mlist = NULL, check_data = TRUE)
```

## Arguments

- ...:

  One or more `brmsfit` objects.

- mlist:

  Optional list of one or more `brmsfit` objects.

- check_data:

  Logical; indicates if the data should be checked for being the same
  across models (defaults to `TRUE`). Setting it to `FALSE` may be
  useful for instance when combining models fitted on multiple imputed
  data sets.

## Value

A `brmsfit` object.

## Details

This function just takes the first model and replaces its `stanfit`
object (slot `fit`) by the combined `stanfit` objects of all models.
