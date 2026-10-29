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
#> Chain 1: Gradient evaluation took 0.000235 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.35 seconds.
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
#> Chain 1:  Elapsed Time: 2.301 seconds (Warm-up)
#> Chain 1:                2.418 seconds (Sampling)
#> Chain 1:                4.719 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000227 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.27 seconds.
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
#> Chain 2:  Elapsed Time: 2.341 seconds (Warm-up)
#> Chain 2:                2.089 seconds (Sampling)
#> Chain 2:                4.43 seconds (Total)
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
#> Intercept[1]    -0.02      0.06    -0.12     0.09 1.00     2340     1525
#> Intercept[2]     0.86      0.10     0.66     1.06 1.00     2059     1682
#> Intercept[3]    -0.16      0.41    -1.08     0.56 1.00     1264     1222
#> period           0.12      0.09    -0.05     0.30 1.00     2120     1503
#> carry           -0.08      0.09    -0.26     0.10 1.00     1707     1448
#> treat[1]        -0.54      0.14    -0.83    -0.26 1.00     1817     1380
#> treat[2]        -0.36      0.22    -0.79     0.06 1.00     1743     1434
#> treat[3]         0.79      0.82    -0.59     2.55 1.00     1340     1199
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
