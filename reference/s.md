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
#> Chain 1: Gradient evaluation took 6.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.64 seconds.
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
#> Chain 1:  Elapsed Time: 2.758 seconds (Warm-up)
#> Chain 1:                2.199 seconds (Sampling)
#> Chain 1:                4.957 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
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
#> Chain 2:  Elapsed Time: 2.732 seconds (Warm-up)
#> Chain 2:                2.264 seconds (Sampling)
#> Chain 2:                4.996 seconds (Total)
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
#> sds(sx0_1)     3.64      1.78     1.44     8.18 1.00     1346     1483
#> sds(sx1_1)     2.14      1.49     0.27     6.03 1.00      885      978
#> sds(sx2_1)    21.24      6.12    12.04    35.76 1.00      518      908
#> sds(sx3_1)     2.77      2.19     0.18     8.76 1.00      773      811
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     7.67      0.15     7.37     7.96 1.01     3142     1360
#> sx0_1         3.73      6.37    -7.62    17.92 1.00     1174     1391
#> sx1_1        12.34      4.48     4.59    22.02 1.00     1122     1267
#> sx2_1        28.13     15.97    -2.66    60.30 1.00     1168     1210
#> sx3_1         7.00      7.27    -2.79    24.73 1.01      889      815
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.09      0.12     1.88     2.34 1.00     2600     1149
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
#> Chain 1: Gradient evaluation took 3.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 1:  Elapsed Time: 3.389 seconds (Warm-up)
#> Chain 1:                4.328 seconds (Sampling)
#> Chain 1:                7.717 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.24 seconds.
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
#> Chain 2:  Elapsed Time: 3.564 seconds (Warm-up)
#> Chain 2:                4.382 seconds (Sampling)
#> Chain 2:                7.946 seconds (Total)
#> Chain 2: 
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit2)
#> Warning: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ t2(x0, x1) + s(x2, by = x3) 
#>    Data: dat (Number of observations: 200) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Smoothing Spline Hyperparameters:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sds(t2x0x1_1)     2.58      2.16     0.11     7.88 1.00     1296     1326
#> sds(t2x0x1_2)     3.92      2.77     0.24    10.77 1.00      823      688
#> sds(t2x0x1_3)     6.92      3.33     1.55    14.55 1.00      979      721
#> sds(sx2x3_1)     32.11     11.44    15.64    59.04 1.00      588     1027
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     6.58      0.39     5.80     7.34 1.00     2909     1671
#> t2x0x1_1      1.93      0.24     1.46     2.40 1.00     2242     1559
#> t2x0x1_2     -0.42      0.24    -0.90     0.04 1.00     2355     1563
#> t2x0x1_3      0.19      0.27    -0.32     0.73 1.00     2771     1887
#> sx2:x3_1     25.96     30.45   -35.46    82.13 1.00      826     1124
#> sx2:x3_2     24.68     18.49   -11.78    61.11 1.00     1457     1251
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.65      0.14     2.39     2.95 1.00     2821     1395
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_smooths(fit2), ask = FALSE)


# }
```
