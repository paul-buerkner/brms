# Draws of a Distributional Parameter

Get draws of a distributional parameter from a `brmsprep` or
`mvbrmsprep` object. This function is primarily useful when developing
custom families or packages depending on brms. This function lets
callers easily handle both the case when the distributional parameter is
predicted directly, via a (non-)linear predictor or fixed to a constant.
See the vignette
[`vignette("brms_customfamilies")`](https://paulbuerkner.com/brms/articles/brms_customfamilies.md)
for an example use case.

## Usage

``` r
get_dpar(prep, dpar, i = NULL, inv_link = NULL)
```

## Arguments

- prep:

  A 'brmsprep' or 'mvbrmsprep' object created by
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md).

- dpar:

  Name of the distributional parameter.

- i:

  The observation numbers for which predictions shall be extracted. If
  `NULL` (the default), all observation will be extracted. Ignored if
  `dpar` is not predicted.

- inv_link:

  Should the inverse link function be applied? If `NULL` (the default),
  the value is chosen internally. In particular, `inv_link` is `TRUE` by
  default for custom families.

## Value

If the parameter is predicted and `i` is `NULL` or `length(i) > 1`, an
`S x N` matrix. If the parameter it not predicted or `length(i) == 1`, a
vector of length `S`. Here `S` is the number of draws and `N` is the
number of observations or length of `i` if specified.

## Examples

``` r
# \dontrun{
posterior_predict_my_dist <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  mypar <- brms::get_dpar(prep, "mypar", i = i)
  my_rng(mu, mypar)
}
# }
```
