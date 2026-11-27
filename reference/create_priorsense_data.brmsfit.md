# Prior sensitivity: Create priorsense data

The `create_priorsense_data.brmsfit` method can be used to create the
data structure needed by the priorsense package for performing
power-scaling sensitivity analysis. This method is called automatically
when performing powerscaling via
[`powerscale`](https://mc-stan.org/priorsense/reference/powerscale-overview.html)
or other related functions, so you will rarely need to call it manually
yourself.

## Usage

``` r
create_priorsense_data.brmsfit(x, ...)
```

## Arguments

- x:

  A `brmsfit` object.

- ...:

  Additional arguments passed to
  [`log_lik`](https://paulbuerkner.com/brms/reference/log_lik.brmsfit.md),
  for example `newdata`.

## Value

A `priorsense_data` object to be used in conjunction with the priorsense
package.

## Examples

``` r
# \dontrun{
# fit a model with non-uniform priors
fit <- brm(rating ~ treat + period + carry,
           data = inhaler, family = sratio(),
           prior = set_prior("normal(0, 0.5)"))
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000133 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.33 seconds.
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
#> Chain 1:  Elapsed Time: 0.594 seconds (Warm-up)
#> Chain 1:                0.536 seconds (Sampling)
#> Chain 1:                1.13 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000109 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.09 seconds.
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
#> Chain 2:  Elapsed Time: 0.597 seconds (Warm-up)
#> Chain 2:                0.591 seconds (Sampling)
#> Chain 2:                1.188 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 8.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.83 seconds.
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
#> Chain 3:  Elapsed Time: 0.607 seconds (Warm-up)
#> Chain 3:                0.558 seconds (Sampling)
#> Chain 3:                1.165 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 9.4e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.94 seconds.
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
#> Chain 4:  Elapsed Time: 0.593 seconds (Warm-up)
#> Chain 4:                0.574 seconds (Sampling)
#> Chain 4:                1.167 seconds (Total)
#> Chain 4: 
summary(fit)
#>  Family: sratio 
#>   Links: mu = logit 
#> Formula: rating ~ treat + period + carry 
#>    Data: inhaler (Number of observations: 572) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept[1]     0.55      0.09     0.37     0.73 1.00     5792     2925
#> Intercept[2]     2.24      0.22     1.82     2.69 1.00     5365     2941
#> Intercept[3]     1.07      0.44     0.24     1.98 1.00     6355     3040
#> treat           -0.64      0.20    -1.04    -0.25 1.00     4194     2928
#> period           0.20      0.16    -0.11     0.51 1.00     5331     2860
#> carry           -0.26      0.15    -0.56     0.04 1.00     3882     2891
#> 
#> Further Distributional Parameters:
#>      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> disc     1.00      0.00     1.00     1.00   NA       NA       NA
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

# The following code requires the 'priorsense' package to be installed:
library(priorsense)

# perform power-scaling of the prior
powerscale(fit, alpha = 1.5, component = "prior")
#> # A draws_df: 1000 iterations, 4 chains, and 10 variables
#>    b_Intercept[1] b_Intercept[2] b_Intercept[3] b_treat b_period b_carry disc
#> 1            0.47            2.7           1.98   -0.77    0.282  -0.201    1
#> 2            0.54            1.8           0.55   -0.54    0.137  -0.239    1
#> 3            0.58            1.8           0.51   -0.42   -0.062  -0.255    1
#> 4            0.51            2.5           1.25   -0.95    0.358  -0.205    1
#> 5            0.48            1.9           0.98   -0.47    0.098  -0.408    1
#> 6            0.63            2.6           1.24   -0.50    0.278  -0.486    1
#> 7            0.66            1.9           0.70   -0.82    0.325  -0.103    1
#> 8            0.59            2.1           1.19   -0.84    0.201  -0.280    1
#> 9            0.51            2.2           0.73   -0.68    0.148  -0.352    1
#> 10           0.35            1.8           0.59   -0.52    0.070   0.027    1
#>    Intercept[1]
#> 1          0.47
#> 2          0.54
#> 3          0.58
#> 4          0.51
#> 5          0.48
#> 6          0.63
#> 7          0.66
#> 8          0.59
#> 9          0.51
#> 10         0.35
#> # ... with 3990 more draws, and 2 more variables
#> # ... hidden reserved variables {'.log_weight', '.chain', '.iteration', '.draw'}
#> 
#> power-scaling
#>  alpha: 1.5 
#>  scaled component: prior 
#>  selection: 
#>  pareto-k: -0.24 
#>  pareto-k threshold: 0.72 
#>  resampled: FALSE 
#>  transform: identity 

# perform power-scaling sensitivity checks
powerscale_sensitivity(fit)
#> Sensitivity based on cjs_dist
#> Prior selection: all priors
#> Likelihood selection: all data
#> 
#>        variable prior likelihood                           diagnosis
#>  b_Intercept[1] 0.011      0.092                                   -
#>  b_Intercept[2] 0.029      0.087                                   -
#>  b_Intercept[3] 0.027      0.088                                   -
#>         b_treat 0.078      0.116 potential prior-likelihood conflict
#>        b_period 0.025      0.081                                   -
#>         b_carry 0.029      0.092                                   -
#>            disc   NaN        NaN                                <NA>
#>    Intercept[1] 0.010      0.092                                   -
#>    Intercept[2] 0.029      0.087                                   -
#>    Intercept[3] 0.027      0.088                                   -

# create power-scaling sensitivity plots (for one variable)
powerscale_plot_dens(fit, variable = "b_treat")
#> 
#> ── Power-scaling sensitivity density plot: ──
#> 
#> The plot shows posterior density depending on the degree of power-scaling
#> (alpha). Overlapping lines indicate low sensitivity. Wider gaps between lines
#> indicate greater sensitivity. Estimates with high Pareto k (dashed lines) may
#> be inaccurate.
#> Disable this help text with `help_text = FALSE` or
#> `options(priorsense.plot_help_text = FALSE)`

# }
```
