# Transform `brmsfit` to `draws` objects

Transform a `brmsfit` object to a format supported by the posterior
package.

## Usage

``` r
# S3 method for class 'brmsfit'
as_draws(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'brmsfit'
as_draws_matrix(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'brmsfit'
as_draws_array(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'brmsfit'
as_draws_df(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'brmsfit'
as_draws_list(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'brmsfit'
as_draws_rvars(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)
```

## Arguments

- x:

  A `brmsfit` object or another R object for which the methods are
  defined.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

- regex:

  Logical; Should variable should be treated as a (vector of) regular
  expressions? Any variable in `x` matching at least one of the regular
  expressions will be selected. Defaults to `FALSE`.

- inc_warmup:

  Should warmup draws be included? Defaults to `FALSE`.

- ...:

  Arguments passed to individual methods (if applicable).

## Details

To subset iterations, chains, or draws, use the
[`subset_draws`](https://mc-stan.org/posterior/reference/subset_draws.html)
method after transforming the `brmsfit` to a `draws` object.

## See also

[`draws`](https://mc-stan.org/posterior/reference/draws.html)
[`subset_draws`](https://mc-stan.org/posterior/reference/subset_draws.html)

## Examples

``` r
# \dontrun{
fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
           data = epilepsy, family = poisson())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.28 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.922 seconds (Warm-up)
#> Chain 1:                1.392 seconds (Sampling)
#> Chain 1:                3.314 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.3e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 1.839 seconds (Warm-up)
#> Chain 2:                1.352 seconds (Sampling)
#> Chain 2:                3.191 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 2.081 seconds (Warm-up)
#> Chain 3:                1.368 seconds (Sampling)
#> Chain 3:                3.449 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.9e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.29 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 2.02 seconds (Warm-up)
#> Chain 4:                1.4 seconds (Sampling)
#> Chain 4:                3.42 seconds (Total)
#> Chain 4: 

# extract posterior draws in an array format
(draws_fit <- as_draws_array(fit))
#> # A draws_array: 1000 iterations, 4 chains, and 68 variables
#> , , variable = b_Intercept
#> 
#>          chain
#> iteration   1   2   3   4
#>         1 1.7 1.8 1.9 2.0
#>         2 1.6 1.7 1.8 1.9
#>         3 1.6 1.7 1.8 2.0
#>         4 1.6 1.7 1.9 1.9
#>         5 1.7 1.6 1.8 1.7
#> 
#> , , variable = b_zAge
#> 
#>          chain
#> iteration     1      2    3       4
#>         1 0.075  0.068 0.14 -0.0318
#>         2 0.120  0.079 0.21 -0.0426
#>         3 0.137 -0.020 0.15 -0.0056
#>         4 0.121  0.063 0.17 -0.0186
#>         5 0.125  0.049 0.23  0.1551
#> 
#> , , variable = b_zBase
#> 
#>          chain
#> iteration    1    2    3    4
#>         1 0.81 0.74 0.65 0.69
#>         2 0.88 0.79 0.79 0.77
#>         3 0.87 0.72 0.78 0.72
#>         4 0.85 0.70 0.82 0.87
#>         5 0.75 0.51 0.85 0.74
#> 
#> , , variable = b_Trt1
#> 
#>          chain
#> iteration     1     2     3     4
#>         1 -0.43 -0.31 -0.35 -0.41
#>         2 -0.23 -0.13 -0.26 -0.44
#>         3 -0.15 -0.20 -0.33 -0.54
#>         4 -0.32 -0.16 -0.24 -0.43
#>         5 -0.22  0.12 -0.25 -0.27
#> 
#> # ... with 995 more iterations, and 64 more variables
posterior::summarize_draws(draws_fit)
#> # A tibble: 68 × 10
#>    variable            mean  median     sd    mad      q5     q95  rhat ess_bulk
#>    <chr>              <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <dbl>    <dbl>
#>  1 b_Intercept       1.77    1.78   0.114  0.111   1.58    1.96    1.01     689.
#>  2 b_zAge            0.0954  0.0925 0.0860 0.0859 -0.0426  0.239   1.01     819.
#>  3 b_zBase           0.696   0.697  0.118  0.117   0.500   0.888   1.00     730.
#>  4 b_Trt1           -0.269  -0.268  0.162  0.163  -0.538  -0.0141  1.01     744.
#>  5 b_zBase:Trt1      0.0604  0.0625 0.159  0.160  -0.198   0.326   1.00     933.
#>  6 sd_patient__Int…  0.583   0.577  0.0704 0.0691  0.478   0.708   1.00     779.
#>  7 Intercept         1.63    1.63   0.0812 0.0816  1.50    1.76    1.00     615.
#>  8 r_patient[1,Int… -0.0516 -0.0485 0.270  0.275  -0.502   0.379   1.00    2504.
#>  9 r_patient[2,Int… -0.0340 -0.0293 0.267  0.271  -0.480   0.388   1.00    2583.
#> 10 r_patient[3,Int… -0.0713 -0.0732 0.300  0.301  -0.572   0.412   1.00    2793.
#> # ℹ 58 more rows
#> # ℹ 1 more variable: ess_tail <dbl>

# extract only certain variables
as_draws_array(fit, variable = "r_patient")
#> # A draws_array: 1000 iterations, 4 chains, and 59 variables
#> , , variable = r_patient[1,Intercept]
#> 
#>          chain
#> iteration      1      2     3     4
#>         1  0.419  0.077 -0.26 -0.27
#>         2  0.085  0.135  0.20 -0.15
#>         3  0.063  0.501  0.03 -0.22
#>         4 -0.095 -0.105 -0.46  0.11
#>         5 -0.098  0.024 -0.64 -0.16
#> 
#> , , variable = r_patient[2,Intercept]
#> 
#>          chain
#> iteration    1      2      3      4
#>         1 0.18  0.087  0.045 -0.636
#>         2 0.35 -0.059 -0.120  0.354
#>         3 0.41  0.329  0.091  0.167
#>         4 0.27 -0.095  0.249 -0.308
#>         5 0.10  0.098 -0.201  0.032
#> 
#> , , variable = r_patient[3,Intercept]
#> 
#>          chain
#> iteration    1     2      3      4
#>         1 0.29 -0.23 -0.512 -0.087
#>         2 0.15  0.26 -0.088 -0.147
#>         3 0.21 -0.28  0.141 -0.120
#>         4 0.68 -0.03  0.167 -0.058
#>         5 0.65 -0.14  0.126  0.112
#> 
#> , , variable = r_patient[4,Intercept]
#> 
#>          chain
#> iteration      1      2      3       4
#>         1 -0.025 -0.085 -0.289  0.0640
#>         2  0.155  0.124 -0.054 -0.0196
#>         3 -0.159  0.354 -0.329  0.0728
#>         4 -0.119  0.121 -0.033 -0.6475
#>         5  0.221  0.031 -0.118 -0.0069
#> 
#> # ... with 995 more iterations, and 55 more variables
as_draws_array(fit, variable = "^b_", regex = TRUE)
#> # A draws_array: 1000 iterations, 4 chains, and 5 variables
#> , , variable = b_Intercept
#> 
#>          chain
#> iteration   1   2   3   4
#>         1 1.7 1.8 1.9 2.0
#>         2 1.6 1.7 1.8 1.9
#>         3 1.6 1.7 1.8 2.0
#>         4 1.6 1.7 1.9 1.9
#>         5 1.7 1.6 1.8 1.7
#> 
#> , , variable = b_zAge
#> 
#>          chain
#> iteration     1      2    3       4
#>         1 0.075  0.068 0.14 -0.0318
#>         2 0.120  0.079 0.21 -0.0426
#>         3 0.137 -0.020 0.15 -0.0056
#>         4 0.121  0.063 0.17 -0.0186
#>         5 0.125  0.049 0.23  0.1551
#> 
#> , , variable = b_zBase
#> 
#>          chain
#> iteration    1    2    3    4
#>         1 0.81 0.74 0.65 0.69
#>         2 0.88 0.79 0.79 0.77
#>         3 0.87 0.72 0.78 0.72
#>         4 0.85 0.70 0.82 0.87
#>         5 0.75 0.51 0.85 0.74
#> 
#> , , variable = b_Trt1
#> 
#>          chain
#> iteration     1     2     3     4
#>         1 -0.43 -0.31 -0.35 -0.41
#>         2 -0.23 -0.13 -0.26 -0.44
#>         3 -0.15 -0.20 -0.33 -0.54
#>         4 -0.32 -0.16 -0.24 -0.43
#>         5 -0.22  0.12 -0.25 -0.27
#> 
#> # ... with 995 more iterations, and 1 more variables

# extract posterior draws in a random variables format
as_draws_rvars(fit)
#> # A draws_rvars: 1000 iterations, 4 chains, and 10 variables
#> $b_Intercept: rvar<1000,4>[1] mean ± sd:
#> [1] 1.8 ± 0.11 
#> 
#> $b_zAge: rvar<1000,4>[1] mean ± sd:
#> [1] 0.095 ± 0.086 
#> 
#> $b_zBase: rvar<1000,4>[1] mean ± sd:
#> [1] 0.7 ± 0.12 
#> 
#> $b_Trt1: rvar<1000,4>[1] mean ± sd:
#> [1] -0.27 ± 0.16 
#> 
#> $b_zBase:Trt1: rvar<1000,4>[1] mean ± sd:
#> [1] 0.06 ± 0.16 
#> 
#> $sd_patient__Intercept: rvar<1000,4>[1] mean ± sd:
#> [1] 0.58 ± 0.07 
#> 
#> $Intercept: rvar<1000,4>[1] mean ± sd:
#> [1] 1.6 ± 0.081 
#> 
#> $r_patient: rvar<1000,4>[59,1] mean ± sd:
#>       Intercept      
#>  [1,] -0.0516 ± 0.27 
#>  [2,] -0.0340 ± 0.27 
#>  [3,] -0.0713 ± 0.30 
#>  [4,] -0.1086 ± 0.29 
#>  [5,]  0.0315 ± 0.24 
#>  [6,]  0.0063 ± 0.23 
#>  [7,] -0.2011 ± 0.28 
#>  [8,]  0.5961 ± 0.25 
#>  [9,] -0.0025 ± 0.24 
#> [10,]  0.8063 ± 0.22 
#> [11,]  0.3449 ± 0.21 
#> [12,]  0.2190 ± 0.21 
#> [13,]  0.0080 ± 0.27 
#> [14,]  0.1571 ± 0.20 
#> [15,] -0.4845 ± 0.29 
#> [16,] -0.7345 ± 0.25 
#> [17,] -0.7437 ± 0.32 
#> [18,] -0.4594 ± 0.37 
#> [19,] -0.1631 ± 0.26 
#> [20,] -0.0058 ± 0.27 
#> [21,] -0.0476 ± 0.27 
#> [22,]  0.1162 ± 0.28 
#> [23,] -0.2569 ± 0.27 
#> [24,]  0.3284 ± 0.22 
#> [25,]  1.1310 ± 0.17 
#> [26,] -0.7133 ± 0.35 
#> [27,] -0.1478 ± 0.31 
#> [28,]  0.4598 ± 0.20 
#> [29,] -0.2526 ± 0.25 
#> [30,]  0.1567 ± 0.22 
#> [31,] -0.3691 ± 0.33 
#> [32,]  0.1736 ± 0.29 
#> [33,]  0.4581 ± 0.28 
#> [34,] -0.1982 ± 0.28 
#> [35,]  1.3297 ± 0.16 
#> [36,]  0.4055 ± 0.26 
#> [37,] -0.0248 ± 0.31 
#> [38,] -0.5493 ± 0.25 
#> [39,]  0.2523 ± 0.22 
#> [40,] -0.5377 ± 0.36 
#> [41,] -0.5484 ± 0.33 
#> [42,] -0.0830 ± 0.31 
#> [43,]  0.7476 ± 0.19 
#> [44,]  0.2891 ± 0.23 
#> [45,]  0.4281 ± 0.22 
#> [46,] -0.1836 ± 0.33 
#> [47,]  0.4203 ± 0.20 
#> [48,] -0.6977 ± 0.37 
#> [49,] -0.4728 ± 0.49 
#> [50,] -0.1247 ± 0.28 
#> [51,]  0.1095 ± 0.22 
#> [52,] -0.5858 ± 0.29 
#> [53,]  0.7304 ± 0.19 
#> [54,] -0.2793 ± 0.30 
#> [55,]  0.1335 ± 0.26 
#> [56,]  1.2501 ± 0.19 
#> [57,] -0.6005 ± 0.32 
#> [58,] -1.2529 ± 0.42 
#> [59,] -0.1676 ± 0.31 
#> 
#> # ... with 2 more variables
# }
```
