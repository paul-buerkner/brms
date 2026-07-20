# Linear and Non-linear formulas in brms

Helper functions to specify linear and non-linear formulas for use with
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

## Usage

``` r
nlf(formula, ..., flist = NULL, dpar = NULL, resp = NULL, loop = NULL)

lf(
  ...,
  flist = NULL,
  dpar = NULL,
  resp = NULL,
  center = NULL,
  cmc = NULL,
  sparse = NULL,
  decomp = NULL
)

acformula(autocor, resp = NULL)

set_nl(nl = TRUE, dpar = NULL, resp = NULL)

set_rescor(rescor = TRUE)

set_mecor(mecor = TRUE)
```

## Arguments

- formula:

  Non-linear formula for a distributional parameter. The name of the
  distributional parameter can either be specified on the left-hand side
  of `formula` or via argument `dpar`.

- ...:

  Additional `formula` objects to specify predictors of non-linear and
  distributional parameters. Formulas can either be named directly or
  contain names on their left-hand side. Alternatively, it is possible
  to fix parameters to certain values by passing numbers or character
  strings in which case arguments have to be named to provide the
  parameter names. See 'Details' for more information.

- flist:

  Optional list of formulas, which are treated in the same way as
  formulas passed via the `...` argument.

- dpar:

  Optional character string specifying the distributional parameter to
  which the formulas passed via `...` and `flist` belong.

- resp:

  Optional character string specifying the response variable to which
  the formulas passed via `...` and `flist` belong. Only relevant in
  multivariate models.

- loop:

  Logical; Only used in non-linear models. Indicates if the computation
  of the non-linear formula should be done inside (`TRUE`) or outside
  (`FALSE`) a loop over observations. Defaults to `TRUE`.

- center:

  Logical; Indicates if the population-level design matrix should be
  centered, which usually increases sampling efficiency. See the
  'Details' section for more information. Defaults to `TRUE` for
  distributional parameters and to `FALSE` for non-linear parameters.

- cmc:

  Logical; Indicates whether automatic cell-mean coding should be
  enabled when removing the intercept by adding `0` to the right-hand of
  model formulas. Defaults to `TRUE` to mirror the behavior of standard
  R formula parsing.

- sparse:

  Logical; indicates whether the population-level design matrices should
  be treated as sparse (defaults to `FALSE`). For design matrices with
  many zeros, this can considerably reduce required memory. Sampling
  speed is currently not improved or even slightly decreased.

- decomp:

  Optional name of the decomposition used for the population-level
  design matrix. Defaults to `NULL` that is no decomposition. Other
  options currently available are `"QR"` for the QR decomposition that
  helps in fitting models with highly correlated predictors.

- autocor:

  A one sided formula containing autocorrelation terms. All none
  autocorrelation terms in `autocor` will be silently ignored.

- nl:

  Logical; Indicates whether `formula` should be treated as specifying a
  non-linear model. By default, `formula` is treated as an ordinary
  linear model formula.

- rescor:

  Logical; Indicates if residual correlation between the response
  variables should be modeled. Currently this is only possible in
  multivariate `gaussian` and `student` models. Only relevant in
  multivariate models.

- mecor:

  Logical; Indicates if correlations between latent variables defined by
  [`me`](https://paulbuerkner.com/brms/reference/me.md) terms should be
  modeled. Defaults to `TRUE`.

## Value

For `lf` and `nlf` a `list` that can be passed to
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
or added to an existing `brmsformula` or `mvbrmsformula` object. For
`set_nl` and `set_rescor` a logical value that can be added to an
existing `brmsformula` or `mvbrmsformula` object.

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)

## Examples

``` r
# add more formulas to the model
bf(y ~ 1) +
  nlf(sigma ~ a * exp(b * x)) +
  lf(a ~ x, b ~ z + (1|g)) +
  gaussian()
#> y ~ 1 
#> sigma ~ a * exp(b * x)
#> a ~ x
#> b ~ z + (1 | g)

# specify 'nl' later on
bf(y ~ a * inv_logit(x * b)) +
  lf(a + b ~ z) +
  set_nl(TRUE)
#> y ~ a * inv_logit(x * b) 
#> a ~ z
#> b ~ z

# specify a multivariate model
bf(y1 ~ x + (1|g)) +
  bf(y2 ~ z) +
  set_rescor(TRUE)
#> y1 ~ x + (1 | g) 
#> y2 ~ z 

# add autocorrelation terms
bf(y ~ x) + acformula(~ arma(p = 1, q = 1) + car(W))
#> y ~ x 
#> autocor ~ arma(p = 1, q = 1) + car(W)
```
