# Descriptions of `brmshypothesis` Objects

A `brmshypothesis` object contains posterior draws as well as summary
statistics of non-linear hypotheses as returned by
[`hypothesis`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.md).

## Usage

``` r
# S3 method for class 'brmshypothesis'
print(
  x,
  digits = 2,
  chars = 20,
  short = getOption("brms.short_summary", FALSE),
  ...
)

# S3 method for class 'brmshypothesis'
plot(
  x,
  nvariables = 5,
  N = NULL,
  ignore_prior = FALSE,
  chars = 40,
  colors = NULL,
  theme = NULL,
  ask = TRUE,
  plot = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `brmsfit`.

- digits:

  Minimal number of significant digits, see
  [`print.default`](https://rdrr.io/r/base/print.default.html).

- chars:

  Maximum number of characters of each hypothesis to print or plot. If
  `NULL`, print the full hypotheses. Defaults to `20`.

- short:

  A flag indicating whether to provide a shorter summary with less
  informational text. Defaults to `FALSE`. Can be set globally for the
  current session via the `brms.short_summary` option.

- ...:

  Currently ignored.

- nvariables:

  The number of variables (parameters) plotted per page.

- N:

  Deprecated alias of `nvariables`.

- ignore_prior:

  A flag indicating if prior distributions should also be plotted. Only
  used if priors were specified on the relevant parameters.

- colors:

  Two values specifying the colors of the posterior and prior density
  respectively. If `NULL` (the default) colors are taken from the
  current color scheme of the bayesplot package.

- theme:

  A [`theme`](https://ggplot2.tidyverse.org/reference/theme.html) object
  modifying the appearance of the plots. For some basic themes see
  `ggtheme` and
  [`theme_default`](https://mc-stan.org/bayesplot/reference/theme_default.html).

- ask:

  Logical; indicates if the user is prompted before a new page is
  plotted. Only used if `plot` is `TRUE`.

- plot:

  Logical; indicates if plots should be plotted directly in the active
  graphic device. Defaults to `TRUE`.

## Details

The two most important elements of a `brmshypothesis` object are
`hypothesis`, which is a data.frame containing the summary estimates of
the hypotheses, and `samples`, which is a data.frame containing the
corresponding posterior draws.

## See also

[`hypothesis`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.md)
