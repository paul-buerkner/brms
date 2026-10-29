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
#> Chain 1: Gradient evaluation took 0.00018 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.8 seconds.
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
#> Chain 1:  Elapsed Time: 1.172 seconds (Warm-up)
#> Chain 1:                1.074 seconds (Sampling)
#> Chain 1:                2.246 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000139 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.39 seconds.
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
#> Chain 2:  Elapsed Time: 1.175 seconds (Warm-up)
#> Chain 2:                1.078 seconds (Sampling)
#> Chain 2:                2.253 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.00017 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.7 seconds.
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
#> Chain 3:  Elapsed Time: 1.154 seconds (Warm-up)
#> Chain 3:                1.089 seconds (Sampling)
#> Chain 3:                2.243 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.00014 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 1.4 seconds.
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
#> Chain 4:  Elapsed Time: 1.185 seconds (Warm-up)
#> Chain 4:                1.08 seconds (Sampling)
#> Chain 4:                2.265 seconds (Total)
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
#> Intercept[1]     0.55      0.09     0.37     0.72 1.00     5901     3028
#> Intercept[2]     2.24      0.22     1.82     2.69 1.00     5456     3088
#> Intercept[3]     1.07      0.44     0.23     1.96 1.00     6400     2737
#> treat           -0.65      0.20    -1.04    -0.27 1.00     4236     2828
#> period           0.19      0.16    -0.11     0.51 1.00     4916     2853
#> carry           -0.25      0.15    -0.54     0.03 1.00     4676     3295
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
#> 1            0.48            2.2           1.45   -0.48     0.16  -0.260    1
#> 2            0.51            2.8           0.34   -0.64     0.30  -0.422    1
#> 3            0.47            1.7           1.29   -0.87     0.19  -0.193    1
#> 4            0.53            2.0           1.10   -0.78     0.10  -0.141    1
#> 5            0.47            2.0           1.27   -0.66     0.24  -0.012    1
#> 6            0.49            2.0           1.34   -0.83     0.21  -0.180    1
#> 7            0.60            2.4           0.62   -0.65     0.16  -0.374    1
#> 8            0.73            2.7           1.55   -0.72     0.25  -0.236    1
#> 9            0.77            2.5           2.43   -0.79     0.22  -0.207    1
#> 10           0.58            1.8           0.71   -0.49     0.28  -0.291    1
#>    Intercept[1]
#> 1          0.48
#> 2          0.51
#> 3          0.47
#> 4          0.53
#> 5          0.47
#> 6          0.49
#> 7          0.60
#> 8          0.73
#> 9          0.77
#> 10         0.58
#> # ... with 3990 more draws, and 2 more variables
#> # ... hidden reserved variables {'.log_weight', '.chain', '.iteration', '.draw'}
#> 
#> power-scaling
#>  alpha: 1.5 
#>  scaled component: prior 
#>  selection: 
#>  pareto-k: -0.11 
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
#>  b_Intercept[1] 0.010      0.086                                   -
#>  b_Intercept[2] 0.029      0.091                                   -
#>  b_Intercept[3] 0.027      0.083                                   -
#>         b_treat 0.077      0.117 potential prior-likelihood conflict
#>        b_period 0.024      0.076                                   -
#>         b_carry 0.027      0.085                                   -
#>            disc   NaN        NaN                                <NA>
#>    Intercept[1] 0.010      0.085                                   -
#>    Intercept[2] 0.029      0.091                                   -
#>    Intercept[3] 0.027      0.083                                   -

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
