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
#>   x          y            cond__
#> 1 a -0.2925950 x = a & y = -0.29
#> 2 a  0.3363284  x = a & y = 0.34
#> 3 a  0.9652517  x = a & y = 0.97
#> 4 b -0.2925950 x = b & y = -0.29
#> 5 b  0.3363284  x = b & y = 0.34
#> 6 b  0.9652517  x = b & y = 0.97
```
