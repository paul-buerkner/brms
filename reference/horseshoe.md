# Regularized horseshoe priors in brms

Function used to set up regularized horseshoe priors and related
hierarchical shrinkage priors in brms. The function does not evaluate
its arguments – it exists purely to help set up the model.

## Usage

``` r
horseshoe(
  df = 1,
  scale_global = 1,
  df_global = 1,
  scale_slab = 2,
  df_slab = 4,
  par_ratio = NULL,
  autoscale = TRUE,
  main = FALSE
)
```

## Arguments

- df:

  Degrees of freedom of student-t prior of the local shrinkage
  parameters. Defaults to `1`.

- scale_global:

  Scale of the student-t prior of the global shrinkage parameter.
  Defaults to `1`. In linear models, `scale_global` will internally be
  multiplied by the residual standard deviation parameter `sigma`.

- df_global:

  Degrees of freedom of student-t prior of the global shrinkage
  parameter. Defaults to `1`. If `df_global` is greater `1`, the shape
  of the prior will no longer resemble a horseshoe and it may be more
  appropriately called an hierarchical shrinkage prior in this case.

- scale_slab:

  Scale of the Student-t slab. Defaults to `2`. The original
  unregularized horseshoe prior is obtained by setting `scale_slab` to
  infinite, which we can approximate in practice by setting it to a very
  large real value.

- df_slab:

  Degrees of freedom of the student-t slab. Defaults to `4`.

- par_ratio:

  Ratio of the expected number of non-zero coefficients to the expected
  number of zero coefficients. If specified, `scale_global` is ignored
  and internally computed as `par_ratio / sqrt(N)`, where `N` is the
  total number of observations in the data.

- autoscale:

  Logical; indicating whether the horseshoe prior should be scaled using
  the residual standard deviation `sigma` if possible and sensible
  (defaults to `TRUE`). Autoscaling is not applied for distributional
  parameters or when the model does not contain the parameter `sigma`.

- main:

  Logical (defaults to `FALSE`); only relevant if the horseshoe prior
  spans multiple parameter classes. In this case, only arguments given
  in the single instance where `main` is `TRUE` will be used. Arguments
  given in other instances of the prior will be ignored. See the
  Examples section below.

## Value

A character string obtained by
[`match.call()`](https://rdrr.io/r/base/match.call.html) with additional
arguments.

## Details

The horseshoe prior is a special shrinkage prior initially proposed by
Carvalho et al. (2009). It is symmetric around zero with fat tails and
an infinitely large spike at zero. This makes it ideal for sparse models
that have many regression coefficients, although only a minority of them
is non-zero. The horseshoe prior can be applied on all population-level
effects at once (excluding the intercept) by using
`set_prior("horseshoe(1)")`. The `1` implies that the student-t prior of
the local shrinkage parameters has 1 degrees of freedom. This may,
however, lead to an increased number of divergent transition in Stan.
Accordingly, increasing the degrees of freedom to slightly higher values
(e.g., `3`) may often be a better option, although the prior no longer
resembles a horseshoe in this case. Further, the scale of the global
shrinkage parameter plays an important role in amount of shrinkage
applied. It defaults to `1`, but this may result in too few shrinkage
(Piironen & Vehtari, 2016). It is thus possible to change the scale
using argument `scale_global` of the horseshoe prior, for instance
`horseshoe(1, scale_global = 0.5)`. In linear models, `scale_global`
will internally be multiplied by the residual standard deviation
parameter `sigma`. See Piironen and Vehtari (2016) for recommendations
how to properly set the global scale. The degrees of freedom of the
global shrinkage prior may also be adjusted via argument `df_global`.
Piironen and Vehtari (2017) recommend to specifying the ratio of the
expected number of non-zero coefficients to the expected number of zero
coefficients `par_ratio` rather than `scale_global` directly. As
proposed by Piironen and Vehtari (2017), an additional regularization is
applied that only affects non-zero coefficients. The amount of
regularization can be controlled via `scale_slab` and `df_slab`. To make
sure that shrinkage can equally affect all coefficients, predictors
should be one the same scale. Generally, models with horseshoe priors a
more likely than other models to have divergent transitions so that
increasing `adapt_delta` from `0.8` to values closer to `1` will often
be necessary. See the documentation of
[`brm`](https://paulbuerkner.com/brms/reference/brm.md) for instructions
on how to increase `adapt_delta`.

The prior does not account for scale differences of the terms it is
applied on. Accordingly, please make sure that all these terms have a
comparable scale to ensure that shrinkage is applied properly.

Currently, the following classes support the horseshoe prior: `b`
(overall regression coefficients), `sds` (SDs of smoothing splines),
`sdgp` (SDs of Gaussian processes), `ar` (autoregressive coefficients),
`ma` (moving average coefficients), `sderr` (SD of latent residuals),
`sdcar` (SD of spatial CAR structures), `sd` (SD of varying
coefficients).

## References

Carvalho, C. M., Polson, N. G., & Scott, J. G. (2009). Handling sparsity
via the horseshoe. Artificial Intelligence and Statistics.
<http://proceedings.mlr.press/v5/carvalho09a>

Piironen J. & Vehtari A. (2017). On the Hyperprior Choice for the Global
Shrinkage Parameter in the Horseshoe Prior. Artificial Intelligence and
Statistics. <https://arxiv.org/pdf/1610.05559v1>

Piironen, J., and Vehtari, A. (2017). Sparsity information and
regularization in the horseshoe and other shrinkage priors. Electronic
Journal of Statistics. <https://arxiv.org/abs/1707.01694>

## See also

[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md)

## Examples

``` r
set_prior(horseshoe(df = 3, par_ratio = 0.1))
#> b ~ horseshoe(df = 3, par_ratio = 0.1)

# specify the horseshoe prior across multiple parameter classes
set_prior(horseshoe(df = 3, par_ratio = 0.1, main = TRUE), class = "b") +
  set_prior(horseshoe(), class = "sd")
#>                                            prior class coef group resp dpar
#>  horseshoe(df = 3, par_ratio = 0.1, main = TRUE)     b                     
#>                                      horseshoe()    sd                     
#>  nlpar   lb   ub tag source
#>        <NA> <NA>       user
#>        <NA> <NA>       user
```
