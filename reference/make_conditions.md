# Prepare Fully Crossed Conditions

This is a helper function to prepare fully crossed conditions primarily
for use with the `conditions` argument of
[`conditional_effects`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.md).
Automatically creates labels for each row in the `cond__` column.

## Usage

``` r
make_conditions(x, vars, ...)
```

## Arguments

- x:

  An R object from which to extract the variables that should be part of
  the conditions.

- vars:

  Names of the variables that should be part of the conditions.

- ...:

  Arguments passed to
  [`rows2labels`](https://paulbuerkner.com/brms/reference/rows2labels.md).

## Value

A `data.frame` where each row indicates a condition.

## Details

For factor like variables, all levels are used as conditions. For
numeric variables, `mean + (-1:1) * SD` are used as conditions.

## See also

[`conditional_effects`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.md),
[`rows2labels`](https://paulbuerkner.com/brms/reference/rows2labels.md)

## Examples

``` r
df <- data.frame(x = c("a", "b"), y = rnorm(10))
make_conditions(df, vars = c("x", "y"))
#>   x         y            cond__
#> 1 a -0.772919 x = a & y = -0.77
#> 2 a -0.023671 x = a & y = -0.02
#> 3 a  0.725577  x = a & y = 0.73
#> 4 b -0.772919 x = b & y = -0.77
#> 5 b -0.023671 x = b & y = -0.02
#> 6 b  0.725577  x = b & y = 0.73
```
