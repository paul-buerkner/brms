# Category Specific Predictors in brms Models

Category Specific Predictors in brms Models

## Usage

``` r
cs(expr)
```

## Arguments

- expr:

  Expression containing predictors, for which category specific effects
  should be estimated. For evaluation, R formula syntax is applied.

## Details

For detailed documentation see
[`help(brmsformula)`](https://paulbuerkner.com/brms/reference/brmsformula.md)
as well as `vignette("brms_overview")`.

This function is almost solely useful when called in formulas passed to
the brms package.

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)

## Examples

``` r
# \dontrun{
fit <- brm(rating ~ period + carry + cs(treat),
           data = inhaler, family = sratio("cloglog"),
           prior = set_prior("normal(0,5)"), chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000263 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.63 seconds.
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
#> Chain 1:  Elapsed Time: 2.395 seconds (Warm-up)
#> Chain 1:                2.295 seconds (Sampling)
#> Chain 1:                4.69 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000207 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.07 seconds.
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
#> Chain 2:  Elapsed Time: 2.379 seconds (Warm-up)
#> Chain 2:                2.04 seconds (Sampling)
#> Chain 2:                4.419 seconds (Total)
#> Chain 2: 
summary(fit)
#>  Family: sratio 
#>   Links: mu = cloglog 
#> Formula: rating ~ period + carry + cs(treat) 
#>    Data: inhaler (Number of observations: 572) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept[1]    -0.02      0.06    -0.13     0.08 1.00     1956     1394
#> Intercept[2]     0.86      0.10     0.68     1.05 1.00     1965     1419
#> Intercept[3]    -0.18      0.41    -1.06     0.51 1.00     1416     1124
#> period           0.12      0.09    -0.06     0.30 1.00     2130     1603
#> carry           -0.09      0.09    -0.27     0.10 1.00     1462     1466
#> treat[1]        -0.53      0.15    -0.83    -0.24 1.00     1449     1369
#> treat[2]        -0.37      0.21    -0.78     0.03 1.00     1452     1543
#> treat[3]         0.83      0.83    -0.62     2.62 1.00     1492     1055
#> 
#> Further Distributional Parameters:
#>      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> disc     1.00      0.00     1.00     1.00   NA       NA       NA
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(fit, ask = FALSE)


# }
```
