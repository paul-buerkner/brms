# Print a summary for a fitted model represented by a `brmsfit` object

Print a summary for a fitted model represented by a `brmsfit` object

## Usage

``` r
# S3 method for class 'brmsfit'
print(x, digits = 2, short = getOption("brms.short_summary", FALSE), ...)
```

## Arguments

- x:

  An object of class `brmsfit`

- digits:

  The number of significant digits for printing out the summary;
  defaults to 2. The effective sample size is always rounded to
  integers.

- short:

  A flag indicating whether to provide a shorter summary with less
  informational text. Defaults to `FALSE`. Can be set globally for the
  current session via the `brms.short_summary` option.

- ...:

  Additional arguments that would be passed to method `summary` of
  `brmsfit`.

## See also

[`summary.brmsfit`](https://paulbuerkner.com/brms/reference/summary.brmsfit.md)
