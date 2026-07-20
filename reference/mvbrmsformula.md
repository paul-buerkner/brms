# Set up a multivariate model formula for use in brms

Set up a multivariate model formula for use in the brms package allowing
to define (potentially non-linear) additive multilevel models for all
parameters of the assumed response distributions.

## Usage

``` r
mvbrmsformula(..., flist = NULL, rescor = NULL)
```

## Arguments

- ...:

  Objects of class `formula` or `brmsformula`, each specifying a
  univariate model. See
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  for details on how to specify univariate models.

- flist:

  Optional list of formulas, which are treated in the same way as
  formulas passed via the `...` argument.

- rescor:

  Logical; Indicates if residual correlation between the response
  variables should be modeled. Currently, this is only possible in
  multivariate `gaussian` and `student` models. If `NULL` (the default),
  `rescor` is internally set to `TRUE` when possible.

## Value

An object of class `mvbrmsformula`, which is essentially a `list`
containing all model formulas as well as some additional information for
multivariate models.

## Details

See
[`vignette("brms_multivariate")`](https://paulbuerkner.com/brms/articles/brms_multivariate.md)
for a case study.

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`brmsformula-helpers`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.md)

## Examples

``` r
bf1 <- bf(y1 ~ x + (1|g))
bf2 <- bf(y2 ~ s(z))
mvbf(bf1, bf2)
#> y1 ~ x + (1 | g) 
#> y2 ~ s(z) 
```
