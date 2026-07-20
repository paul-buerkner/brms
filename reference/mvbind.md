# Bind response variables in multivariate models

Can be used to specify a multivariate brms model within a single
formula. Outside of
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
it just behaves like
[`cbind`](https://amices.org/mice/reference/cbind.html).

## Usage

``` r
mvbind(...)
```

## Arguments

- ...:

  Same as in [`cbind`](https://amices.org/mice/reference/cbind.html)

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)

## Examples

``` r
bf(mvbind(y1, y2) ~ x)
#> y1 ~ x 
#> y2 ~ x 
```
