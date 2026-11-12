# Defining smooths in brms formulas

Functions used in definition of smooth terms within a model formulas.
The function does not evaluate a (spline) smooth - it exists purely to
help set up a model using spline based smooths.

## Usage

``` r
s(...)

t2(...)
```

## Arguments

- ...:

  Arguments passed to [`mgcv::s`](https://rdrr.io/pkg/mgcv/man/s.html)
  or [`mgcv::t2`](https://rdrr.io/pkg/mgcv/man/t2.html).

## Details

The function defined here are just simple wrappers of the respective
functions of the mgcv package. When using them, please cite the
appropriate references obtained via `citation("mgcv")`.

brms uses the "random effects" parameterization of smoothing splines as
explained in [`mgcv::gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html). A
nice tutorial on this topic can be found in Pedersen et al. (2019). The
answers provided in this [Stan discourse
post](https://discourse.mc-stan.org/t/better-priors-non-flat-for-gams-brms/23012/4)
may also be helpful.

## References

Pedersen, E. J., Miller, D. L., Simpson, G. L., & Ross, N. (2019).
Hierarchical generalized additive models in ecology: an introduction
with mgcv. PeerJ.

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`mgcv::s`](https://rdrr.io/pkg/mgcv/man/s.html),
[`mgcv::t2`](https://rdrr.io/pkg/mgcv/man/t2.html)

## Examples

``` r
# \dontrun{
# simulate some data
dat <- mgcv::gamSim(1, n = 200, scale = 2)
#> Gu & Wahba 4 term additive model

# fit univariate smooths for all predictors
fit1 <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3),
            data = dat, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.83 seconds.
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
#> Chain 1:  Elapsed Time: 4.426 seconds (Warm-up)
#> Chain 1:                3.14 seconds (Sampling)
#> Chain 1:                7.566 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.4 seconds.
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
#> Chain 2:  Elapsed Time: 4.332 seconds (Warm-up)
#> Chain 2:                3.14 seconds (Sampling)
#> Chain 2:                7.472 seconds (Total)
#> Chain 2: 
#> Warning: There were 1 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit1)
#> Warning: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ s(x0) + s(x1) + s(x2) + s(x3) 
#>    Data: dat (Number of observations: 200) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Smoothing Spline Hyperparameters:
#>            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sds(sx0_1)     2.65      1.56     0.85     6.12 1.00     1368     1232
#> sds(sx1_1)     1.84      1.49     0.08     5.62 1.00      918      986
#> sds(sx2_1)    21.34      6.18    12.26    35.73 1.00      654      876
#> sds(sx3_1)     1.78      1.53     0.07     5.55 1.00      683     1036
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     8.11      0.14     7.85     8.37 1.00     3294     1106
#> sx0_1        -0.27      5.13   -11.26     9.86 1.00     1132     1067
#> sx1_1        15.35      4.34     6.77    25.64 1.00     1039      909
#> sx2_1         8.75     15.35   -23.01    36.84 1.00     1017     1336
#> sx3_1        -3.35      4.91   -16.02     3.71 1.00     1005     1045
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.00      0.11     1.79     2.24 1.00     3668     1587
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_smooths(fit1), ask = FALSE)





# fit a more complicated smooth model
fit2 <- brm(y ~ t2(x0, x1) + s(x2, by = x3),
            data = dat, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.61 seconds.
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
#> Chain 1:  Elapsed Time: 4.534 seconds (Warm-up)
#> Chain 1:                5.751 seconds (Sampling)
#> Chain 1:                10.285 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.41 seconds.
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
#> Chain 2:  Elapsed Time: 4.817 seconds (Warm-up)
#> Chain 2:                3.146 seconds (Sampling)
#> Chain 2:                7.963 seconds (Total)
#> Chain 2: 
#> Warning: There were 20 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit2)
#> Warning: There were 20 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ t2(x0, x1) + s(x2, by = x3) 
#>    Data: dat (Number of observations: 200) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Smoothing Spline Hyperparameters:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sds(t2x0x1_1)     3.26      2.27     0.20     8.32 1.01      875      712
#> sds(t2x0x1_2)     2.11      1.86     0.06     6.79 1.00     1212      534
#> sds(t2x0x1_3)     6.70      3.18     1.98    14.29 1.00      752      617
#> sds(sx2x3_1)     24.93      8.20    13.24    45.85 1.00      437      538
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     8.08      0.36     7.35     8.81 1.00     2142     1582
#> t2x0x1_1      2.06      0.21     1.64     2.46 1.00     1811     1201
#> t2x0x1_2      0.03      0.22    -0.40     0.45 1.00     2546     1616
#> t2x0x1_3     -0.00      0.25    -0.48     0.49 1.00     2120     1535
#> sx2:x3_1     28.30     18.01    -9.54    61.98 1.00      949     1145
#> sx2:x3_2     10.30     13.32   -17.75    34.50 1.00     1055     1171
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.54      0.14     2.28     2.82 1.00     2337     1619
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_smooths(fit2), ask = FALSE)


# }
```
