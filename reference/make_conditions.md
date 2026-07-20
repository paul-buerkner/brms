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
#>   x           y            cond__
#> 1 a -0.83573986 x = a & y = -0.84
#> 2 a  0.05411511  x = a & y = 0.05
#> 3 a  0.94397007  x = a & y = 0.94
#> 4 b -0.83573986 x = b & y = -0.84
#> 5 b  0.05411511  x = b & y = 0.05
#> 6 b  0.94397007  x = b & y = 0.94
```
