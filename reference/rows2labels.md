# Convert Rows to Labels

Convert information in rows to labels for each row.

## Usage

``` r
rows2labels(x, digits = 2, sep = " & ", incl_vars = TRUE, ...)
```

## Arguments

- x:

  A `data.frame` for which to extract labels.

- digits:

  Minimal number of decimal places shown in the labels of numeric
  variables.

- sep:

  A single character string defining the separator between variables
  used in the labels.

- incl_vars:

  Indicates if variable names should be part of the labels. Defaults to
  `TRUE`.

- ...:

  Currently unused.

## Value

A character vector of the same length as the number of rows of `x`.

## See also

[`make_conditions`](https://paulbuerkner.com/brms/reference/make_conditions.md),
[`conditional_effects`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.md)
