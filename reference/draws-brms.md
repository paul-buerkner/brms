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
#> Chain 1: Gradient evaluation took 3.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.37 seconds.
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
#> Chain 1:  Elapsed Time: 1.982 seconds (Warm-up)
#> Chain 1:                1.553 seconds (Sampling)
#> Chain 1:                3.535 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.8e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.28 seconds.
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
#> Chain 2:  Elapsed Time: 2.096 seconds (Warm-up)
#> Chain 2:                1.385 seconds (Sampling)
#> Chain 2:                3.481 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.6e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
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
#> Chain 3:  Elapsed Time: 1.959 seconds (Warm-up)
#> Chain 3:                1.658 seconds (Sampling)
#> Chain 3:                3.617 seconds (Total)
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
#> Chain 4:  Elapsed Time: 2.057 seconds (Warm-up)
#> Chain 4:                1.653 seconds (Sampling)
#> Chain 4:                3.71 seconds (Total)
#> Chain 4: 

# extract posterior draws in an array format
(draws_fit <- as_draws_array(fit))
#> # A draws_array: 1000 iterations, 4 chains, and 68 variables
#> , , variable = b_Intercept
#> 
#>          chain
#> iteration   1   2   3   4
#>         1 2.0 1.7 1.6 1.7
#>         2 1.9 1.7 1.7 1.7
#>         3 2.0 1.7 1.6 1.7
#>         4 1.4 1.8 1.6 1.7
#>         5 1.6 1.8 1.7 1.7
#> 
#> , , variable = b_zAge
#> 
#>          chain
#> iteration     1       2     3    4
#>         1 0.057 -0.0191 0.183 0.30
#>         2 0.040 -0.0161 0.225 0.37
#>         3 0.017  0.0035 0.234 0.27
#>         4 0.236  0.0521 0.149 0.22
#>         5 0.237  0.0299 0.075 0.21
#> 
#> , , variable = b_zBase
#> 
#>          chain
#> iteration    1    2    3    4
#>         1 0.67 0.57 0.64 0.67
#>         2 0.79 0.58 0.69 0.70
#>         3 0.71 0.66 0.78 0.60
#>         4 0.63 0.68 0.70 0.62
#>         5 0.66 0.67 0.79 0.65
#> 
#> , , variable = b_Trt1
#> 
#>          chain
#> iteration     1     2      3     4
#>         1 -0.54 -0.14  0.023 -0.23
#>         2 -0.62 -0.23 -0.026 -0.32
#>         3 -0.69 -0.10 -0.039 -0.33
#>         4  0.25 -0.35 -0.040 -0.34
#>         5  0.20 -0.56 -0.444 -0.29
#> 
#> # ... with 995 more iterations, and 64 more variables
posterior::summarize_draws(draws_fit)
#> # A tibble: 68 × 10
#>    variable    mean  median     sd    mad      q5    q95  rhat ess_bulk ess_tail
#>    <chr>      <dbl>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#>  1 b_Inter…  1.76    1.76   0.120  0.116   1.56   1.95    1.00     670.    1117.
#>  2 b_zAge    0.0863  0.0879 0.0900 0.0878 -0.0646 0.234   1.00     718.    1345.
#>  3 b_zBase   0.699   0.698  0.123  0.118   0.494  0.905   1.00     724.    1343.
#>  4 b_Trt1   -0.256  -0.253  0.167  0.165  -0.533  0.0169  1.00     747.    1091.
#>  5 b_zBase…  0.0569  0.0537 0.168  0.166  -0.214  0.335   1.00     916.    1551.
#>  6 sd_pati…  0.583   0.577  0.0714 0.0721  0.477  0.710   1.01     715.    1523.
#>  7 Interce…  1.62    1.62   0.0825 0.0803  1.49   1.76    1.00     675.    1254.
#>  8 r_patie… -0.0286 -0.0254 0.273  0.264  -0.481  0.419   1.00    1993.    2433.
#>  9 r_patie… -0.0216 -0.0189 0.277  0.280  -0.481  0.427   1.00    1594.    2398.
#> 10 r_patie… -0.0560 -0.0491 0.301  0.302  -0.564  0.427   1.00    2211.    2340.
#> # ℹ 58 more rows

# extract only certain variables
as_draws_array(fit, variable = "r_patient")
#> # A draws_array: 1000 iterations, 4 chains, and 59 variables
#> , , variable = r_patient[1,Intercept]
#> 
#>          chain
#> iteration      1      2      3      4
#>         1 -0.062  0.079 -0.048 -0.212
#>         2 -0.219  0.059 -0.018 -0.374
#>         3 -0.125  0.249 -0.211 -0.236
#>         4 -0.018 -0.276  0.484 -0.086
#>         5 -0.084  0.298 -0.437 -0.111
#> 
#> , , variable = r_patient[2,Intercept]
#> 
#>          chain
#> iteration      1      2       3      4
#>         1 -0.224  0.201 -0.0422  0.047
#>         2 -0.083  0.034 -0.0054  0.091
#>         3 -0.215  0.316  0.2369 -0.086
#>         4  0.114 -0.494 -0.0863 -0.411
#>         5  0.020  0.128  0.2596  0.420
#> 
#> , , variable = r_patient[3,Intercept]
#> 
#>          chain
#> iteration      1      2     3      4
#>         1  0.149 -0.531  0.05  0.241
#>         2 -0.488 -0.442  0.26  0.018
#>         3  0.021 -0.325  0.29  0.113
#>         4  0.640  0.190  0.22 -0.058
#>         5 -0.232 -0.067 -0.10 -0.399
#> 
#> , , variable = r_patient[4,Intercept]
#> 
#>          chain
#> iteration      1      2      3      4
#>         1 -0.076  0.291  0.295 -0.530
#>         2 -0.017  0.011  0.496 -0.099
#>         3 -0.269 -0.100 -0.248 -0.618
#>         4 -0.072  0.081 -0.047  0.342
#>         5  0.178 -0.263 -0.032 -0.320
#> 
#> # ... with 995 more iterations, and 55 more variables
as_draws_array(fit, variable = "^b_", regex = TRUE)
#> # A draws_array: 1000 iterations, 4 chains, and 5 variables
#> , , variable = b_Intercept
#> 
#>          chain
#> iteration   1   2   3   4
#>         1 2.0 1.7 1.6 1.7
#>         2 1.9 1.7 1.7 1.7
#>         3 2.0 1.7 1.6 1.7
#>         4 1.4 1.8 1.6 1.7
#>         5 1.6 1.8 1.7 1.7
#> 
#> , , variable = b_zAge
#> 
#>          chain
#> iteration     1       2     3    4
#>         1 0.057 -0.0191 0.183 0.30
#>         2 0.040 -0.0161 0.225 0.37
#>         3 0.017  0.0035 0.234 0.27
#>         4 0.236  0.0521 0.149 0.22
#>         5 0.237  0.0299 0.075 0.21
#> 
#> , , variable = b_zBase
#> 
#>          chain
#> iteration    1    2    3    4
#>         1 0.67 0.57 0.64 0.67
#>         2 0.79 0.58 0.69 0.70
#>         3 0.71 0.66 0.78 0.60
#>         4 0.63 0.68 0.70 0.62
#>         5 0.66 0.67 0.79 0.65
#> 
#> , , variable = b_Trt1
#> 
#>          chain
#> iteration     1     2      3     4
#>         1 -0.54 -0.14  0.023 -0.23
#>         2 -0.62 -0.23 -0.026 -0.32
#>         3 -0.69 -0.10 -0.039 -0.33
#>         4  0.25 -0.35 -0.040 -0.34
#>         5  0.20 -0.56 -0.444 -0.29
#> 
#> # ... with 995 more iterations, and 1 more variables

# extract posterior draws in a random variables format
as_draws_rvars(fit)
#> # A draws_rvars: 1000 iterations, 4 chains, and 10 variables
#> $b_Intercept: rvar<1000,4>[1] mean ± sd:
#> [1] 1.8 ± 0.12 
#> 
#> $b_zAge: rvar<1000,4>[1] mean ± sd:
#> [1] 0.086 ± 0.09 
#> 
#> $b_zBase: rvar<1000,4>[1] mean ± sd:
#> [1] 0.7 ± 0.12 
#> 
#> $b_Trt1: rvar<1000,4>[1] mean ± sd:
#> [1] -0.26 ± 0.17 
#> 
#> $b_zBase:Trt1: rvar<1000,4>[1] mean ± sd:
#> [1] 0.057 ± 0.17 
#> 
#> $sd_patient__Intercept: rvar<1000,4>[1] mean ± sd:
#> [1] 0.58 ± 0.071 
#> 
#> $Intercept: rvar<1000,4>[1] mean ± sd:
#> [1] 1.6 ± 0.082 
#> 
#> $r_patient: rvar<1000,4>[59,1] mean ± sd:
#>       Intercept      
#>  [1,] -0.0286 ± 0.27 
#>  [2,] -0.0216 ± 0.28 
#>  [3,] -0.0560 ± 0.30 
#>  [4,] -0.0837 ± 0.29 
#>  [5,]  0.0320 ± 0.25 
#>  [6,]  0.0345 ± 0.22 
#>  [7,] -0.1807 ± 0.28 
#>  [8,]  0.6319 ± 0.25 
#>  [9,]  0.0271 ± 0.25 
#> [10,]  0.8194 ± 0.22 
#> [11,]  0.3731 ± 0.21 
#> [12,]  0.2327 ± 0.22 
#> [13,]  0.0188 ± 0.27 
#> [14,]  0.1832 ± 0.21 
#> [15,] -0.4767 ± 0.30 
#> [16,] -0.7210 ± 0.25 
#> [17,] -0.7268 ± 0.32 
#> [18,] -0.4512 ± 0.39 
#> [19,] -0.1368 ± 0.26 
#> [20,]  0.0021 ± 0.28 
#> [21,] -0.0366 ± 0.28 
#> [22,]  0.1208 ± 0.30 
#> [23,] -0.2337 ± 0.27 
#> [24,]  0.3410 ± 0.22 
#> [25,]  1.1462 ± 0.18 
#> [26,] -0.6750 ± 0.34 
#> [27,] -0.1453 ± 0.32 
#> [28,]  0.4661 ± 0.22 
#> [29,] -0.2718 ± 0.26 
#> [30,]  0.1653 ± 0.23 
#> [31,] -0.3814 ± 0.33 
#> [32,]  0.1768 ± 0.29 
#> [33,]  0.4466 ± 0.28 
#> [34,] -0.2027 ± 0.28 
#> [35,]  1.3362 ± 0.17 
#> [36,]  0.4165 ± 0.25 
#> [37,] -0.0223 ± 0.29 
#> [38,] -0.5616 ± 0.25 
#> [39,]  0.2480 ± 0.21 
#> [40,] -0.5267 ± 0.37 
#> [41,] -0.5567 ± 0.33 
#> [42,] -0.0619 ± 0.31 
#> [43,]  0.7577 ± 0.19 
#> [44,]  0.2852 ± 0.23 
#> [45,]  0.4439 ± 0.21 
#> [46,] -0.1876 ± 0.33 
#> [47,]  0.4200 ± 0.20 
#> [48,] -0.6962 ± 0.38 
#> [49,] -0.4771 ± 0.50 
#> [50,] -0.1098 ± 0.27 
#> [51,]  0.1075 ± 0.22 
#> [52,] -0.5693 ± 0.29 
#> [53,]  0.7200 ± 0.19 
#> [54,] -0.2582 ± 0.30 
#> [55,]  0.1429 ± 0.27 
#> [56,]  1.2499 ± 0.19 
#> [57,] -0.5991 ± 0.34 
#> [58,] -1.2457 ± 0.43 
#> [59,] -0.1560 ± 0.31 
#> 
#> # ... with 2 more variables
# }
```
