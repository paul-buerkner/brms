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
#> Chain 1: Gradient evaluation took 7.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.74 seconds.
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
#> Chain 1:  Elapsed Time: 5.491 seconds (Warm-up)
#> Chain 1:                4.279 seconds (Sampling)
#> Chain 1:                9.77 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 2:  Elapsed Time: 5.673 seconds (Warm-up)
#> Chain 2:                4.938 seconds (Sampling)
#> Chain 2:                10.611 seconds (Total)
#> Chain 2: 
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit1)
#> Warning: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ s(x0) + s(x1) + s(x2) + s(x3) 
#>    Data: dat (Number of observations: 200) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Smoothing Spline Hyperparameters:
#>            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sds(sx0_1)     2.72      1.42     0.87     6.07 1.00      971     1201
#> sds(sx1_1)     1.72      1.35     0.07     5.10 1.00      630      729
#> sds(sx2_1)    22.44      6.81    12.34    38.27 1.00      481      752
#> sds(sx3_1)     1.80      1.55     0.07     5.68 1.00      647     1142
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     8.12      0.14     7.82     8.39 1.00     2716     1223
#> sx0_1        -0.16      5.13    -9.91    10.18 1.00      998      853
#> sx1_1        15.46      4.06     8.19    25.51 1.00      841      570
#> sx2_1         7.36     15.39   -24.12    38.80 1.00      791     1056
#> sx3_1        -3.43      4.90   -15.61     3.53 1.00      851      733
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.00      0.10     1.81     2.22 1.00     2046     1582
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
#> Chain 1: Gradient evaluation took 5.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.53 seconds.
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
#> Chain 1:  Elapsed Time: 6.359 seconds (Warm-up)
#> Chain 1:                5.016 seconds (Sampling)
#> Chain 1:                11.375 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.35 seconds.
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
#> Chain 2:  Elapsed Time: 5.903 seconds (Warm-up)
#> Chain 2:                6.973 seconds (Sampling)
#> Chain 2:                12.876 seconds (Total)
#> Chain 2: 
#> Warning: There were 24 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit2)
#> Warning: There were 24 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ t2(x0, x1) + s(x2, by = x3) 
#>    Data: dat (Number of observations: 200) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Smoothing Spline Hyperparameters:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sds(t2x0x1_1)     3.38      2.45     0.19     8.80 1.00      797      945
#> sds(t2x0x1_2)     2.16      1.87     0.11     6.77 1.00     1298      917
#> sds(t2x0x1_3)     6.77      3.20     1.81    13.95 1.00      634      376
#> sds(sx2x3_1)     25.37      8.34    12.72    46.57 1.00      515      501
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     8.07      0.37     7.33     8.81 1.00     1811     1404
#> t2x0x1_1      2.06      0.21     1.65     2.47 1.00     1785     1710
#> t2x0x1_2      0.04      0.22    -0.40     0.45 1.00     1952     1327
#> t2x0x1_3      0.00      0.24    -0.47     0.48 1.00     1975     1314
#> sx2:x3_1     27.86     19.46   -12.54    66.89 1.00      773      499
#> sx2:x3_2     10.43     13.67   -17.86    36.39 1.00      829      572
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.53      0.14     2.26     2.81 1.00     1037      548
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_smooths(fit2), ask = FALSE)


# }
```
