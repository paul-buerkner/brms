# Autocorrelation structures

Specify autocorrelation terms in brms models. Currently supported terms
are [`arma`](https://paulbuerkner.com/brms/reference/arma.md),
[`ar`](https://paulbuerkner.com/brms/reference/ar.md),
[`ma`](https://paulbuerkner.com/brms/reference/ma.md),
[`cosy`](https://paulbuerkner.com/brms/reference/cosy.md),
[`unstr`](https://paulbuerkner.com/brms/reference/unstr.md),
[`sar`](https://paulbuerkner.com/brms/reference/sar.md),
[`car`](https://paulbuerkner.com/brms/reference/car.md), and
[`fcor`](https://paulbuerkner.com/brms/reference/fcor.md). Terms can be
directly specified within the formula, or passed to the `autocor`
argument of
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
in the form of a one-sided formula. For deprecated ways of specifying
autocorrelation terms, see
[`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md).

## Details

The autocor term functions are almost solely useful when called in
formulas passed to the brms package. They do not evaluate its arguments
– but exist purely to help set up a model with autocorrelation terms.

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`acformula`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.md),
[`arma`](https://paulbuerkner.com/brms/reference/arma.md),
[`ar`](https://paulbuerkner.com/brms/reference/ar.md),
[`ma`](https://paulbuerkner.com/brms/reference/ma.md),
[`cosy`](https://paulbuerkner.com/brms/reference/cosy.md),
[`unstr`](https://paulbuerkner.com/brms/reference/unstr.md),
[`sar`](https://paulbuerkner.com/brms/reference/sar.md),
[`car`](https://paulbuerkner.com/brms/reference/car.md),
[`fcor`](https://paulbuerkner.com/brms/reference/fcor.md)

## Examples

``` r
# specify autocor terms within the formula
y ~ x + arma(p = 1, q = 1) + car(M)
#> y ~ x + arma(p = 1, q = 1) + car(M)
#> <environment: 0x55d71f3a8c70>

# specify autocor terms in the 'autocor' argument
bf(y ~ x, autocor = ~ arma(p = 1, q = 1) + car(M))
#> y ~ x 
#> autocor ~ arma(p = 1, q = 1) + car(M)

# specify autocor terms via 'acformula'
bf(y ~ x) + acformula(~ arma(p = 1, q = 1) + car(M))
#> y ~ x 
#> autocor ~ arma(p = 1, q = 1) + car(M)
```
