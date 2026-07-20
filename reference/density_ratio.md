# Compute Density Ratios

Compute the ratio of two densities at given points based on draws of the
corresponding distributions.

## Usage

``` r
density_ratio(x, y = NULL, point = 0, n = 4096, ...)
```

## Arguments

- x:

  Vector of draws from the first distribution, usually the posterior
  distribution of the quantity of interest.

- y:

  Optional vector of draws from the second distribution, usually the
  prior distribution of the quantity of interest. If `NULL` (the
  default), only the density of `x` will be evaluated.

- point:

  Numeric values at which to evaluate and compare the densities.
  Defaults to `0`.

- n:

  Single numeric value. Influences the accuracy of the density
  estimation. See [`density`](https://rdrr.io/r/stats/density.html) for
  details.

- ...:

  Further arguments passed to
  [`density`](https://rdrr.io/r/stats/density.html).

## Value

A vector of length equal to `length(point)`. If `y` is provided, the
density ratio of `x` against `y` is returned. Else, only the density of
`x` is returned.

## Details

In order to achieve sufficient accuracy in the density estimation, more
draws than usual are required. That is you may need an effective sample
size of 10,000 or more to reliably estimate the densities.

## Examples

``` r
x <- rnorm(10000)
y <- rnorm(10000, mean = 1)
density_ratio(x, y, point = c(0, 1))
#> [1] 1.556592 0.637758
```
