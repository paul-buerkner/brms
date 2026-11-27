# Set up Gaussian process terms in brms

Set up a Gaussian process (GP) term in brms. The function does not
evaluate its arguments – it exists purely to help set up a model with GP
terms.

## Usage

``` r
gp(
  ...,
  by = NA,
  k = NA,
  cov = "exp_quad",
  iso = TRUE,
  gr = TRUE,
  cmc = TRUE,
  scale = TRUE,
  c = 5/4
)
```

## Arguments

- ...:

  One or more predictors for the GP.

- by:

  A numeric or factor variable of the same length as each predictor. In
  the numeric vector case, the elements multiply the values returned by
  the GP. In the factor variable case, a separate GP is fitted for each
  factor level.

- k:

  Optional number of basis functions for computing Hilbert-space
  approximate GPs. If `NA` (the default), exact GPs are computed.

- cov:

  Name of the covariance kernel. Currently supported are `"exp_quad"`
  (exponentiated-quadratic kernel; default), `"matern32"` (Matern 3/2
  kernel), `"matern52"` (Matern 5/2 kernel), and `"exponential"`
  (exponential kernel; alias: `"matern12"`).

- iso:

  A flag to indicate whether an isotropic (`TRUE`; the default) or a
  non-isotropic GP should be used. In the former case, the same amount
  of smoothing is applied to all predictors. In the latter case,
  predictors may have different smoothing. Ignored if only a single
  predictor is supplied.

- gr:

  Logical; Indicates if auto-grouping should be used (defaults to
  `TRUE`). If enabled, observations sharing the same predictor values
  will be represented by the same latent variable in the GP. This will
  improve sampling efficiency drastically if the number of unique
  predictor combinations is small relative to the number of
  observations.

- cmc:

  Logical; Only relevant if `by` is a factor. If `TRUE` (the default),
  cell-mean coding is used for the `by`-factor, that is one GP per level
  is estimated. If `FALSE`, contrast GPs are estimated according to the
  contrasts set for the `by`-factor.

- scale:

  Logical; If `TRUE` (the default), predictors are scaled so that the
  maximum Euclidean distance between two points is 1. This often
  improves sampling speed and convergence. Scaling also affects the
  estimated length-scale parameters in that they resemble those of
  scaled predictors (not of the original predictors) if `scale` is
  `TRUE`.

- c:

  Numeric value only used in approximate GPs. Defines the multiplicative
  constant of the predictors' range over which predictions should be
  computed. A good default could be `c = 5/4` but we are still working
  on providing better recommendations.

## Value

An object of class `'gp_term'`, which is a list of arguments to be
interpreted by the formula parsing functions of brms.

## Details

A GP is a stochastic process, which describes the relation between one
or more predictors \\x = (x_1, ..., x_d)\\ and a response \\f(x)\\,
where \\d\\ is the number of predictors. A GP is the generalization of
the multivariate normal distribution to an infinite number of
dimensions. Thus, it can be interpreted as a prior over functions. The
values of \\f( )\\ at any finite set of locations are jointly
multivariate normal, with a covariance matrix defined by the covariance
kernel \\k_p(x_i, x_j)\\, where \\p\\ is the vector of parameters of the
GP: \$\$(f(x_1), \ldots f(x_n) \sim MVN(0, (k_p(x_i, x_j))\_{i,j=1}^n)
.\$\$ The smoothness and general behavior of the function \\f\\ depends
only on the choice of covariance kernel. For a more detailed
introduction to Gaussian processes, see
<https://en.wikipedia.org/wiki/Gaussian_process>.

For mathematical details on the supported kernels, please see the Stan
manual:
<https://mc-stan.org/docs/functions-reference/matrix_operations.html>
under "Gaussian Process Covariance Functions".

There are several parameters estimated for GPs, the most important of
which are as follows:

|  |  |  |  |
|----|----|----|----|
| **Parameter** | **Notation** | **Support** | **Meaning** |
| `lscale` | \\\ell\\ | \\\mathbb{R}^+\\ | length-scale of the GP's covariance kernel |
| `sdgp` | \\\sigma\\ | \\\mathbb{R}^+\\ | marginal standard deviation of the GP's covariance kernel |
| `zgp` | \\z\\ | \\\mathbb{R}\\ | latent variable values of the training data |

Assuming an exponentiated quadratic covariance structure, the parameters
can be broadly interpreted as follows (see also
<https://mc-stan.org/docs/stan-users-guide/gaussian-processes.html#gaussian-process-regression>).
Note that in the above documentation the parameter \\\sigma\\ is denoted
as \\\alpha\\ instead, and the parameter \\\ell\\ as \\\rho\\. The
length-scale \\\ell\\ controls the frequency of the functions
represented by the GP prior i.e., values of \\\ell \gg 0\\ lead to
lower-frequency functions, while values of \\\ell \approx 0\\ lead to
higher-frequency functions. In slightly simpler terms, \\\ell\\ sets the
distance over which observations in the input space are strongly
correlated. The marginal standard deviation \\\sigma\\ controls the
magnitude of the range of the function represented by the GP i.e., it
represents how much the values of the function tend to deviate from the
mean level. Lastly, the parameter \\z\\ represents latent variable
values per observation (or unique input point).

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)

## Examples

``` r
# \dontrun{
# simulate data using the mgcv package
dat <- mgcv::gamSim(1, n = 30, scale = 2)
#> Gu & Wahba 4 term additive model

# fit a simple GP model
fit1 <- brm(y ~ gp(x2), dat, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 1:  Elapsed Time: 1.156 seconds (Warm-up)
#> Chain 1:                0.8 seconds (Sampling)
#> Chain 1:                1.956 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.32 seconds.
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
#> Chain 2:  Elapsed Time: 1.136 seconds (Warm-up)
#> Chain 2:                1.069 seconds (Sampling)
#> Chain 2:                2.205 seconds (Total)
#> Chain 2: 
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit1)
#> Warning: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ gp(x2) 
#>    Data: dat (Number of observations: 30) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Gaussian Process Hyperparameters:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sdgp(gpx2)       3.82      1.51     1.96     7.79 1.01      485      470
#> lscale(gpx2)     0.12      0.07     0.05     0.27 1.01      205      540
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     7.53      1.83     3.61    10.89 1.00      770      806
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.05      0.37     1.48     2.91 1.00      976     1456
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
me1 <- conditional_effects(fit1, ndraws = 200, spaghetti = TRUE)
plot(me1, ask = FALSE, points = TRUE)


# fit a more complicated GP model and use an approximate GP for x2
fit2 <- brm(y ~ gp(x0) + x1 + gp(x2, k = 10) + x3, dat, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.58 seconds.
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
#> Chain 1:  Elapsed Time: 1.341 seconds (Warm-up)
#> Chain 1:                1.769 seconds (Sampling)
#> Chain 1:                3.11 seconds (Total)
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
#> Chain 2:  Elapsed Time: 1.408 seconds (Warm-up)
#> Chain 2:                0.916 seconds (Sampling)
#> Chain 2:                2.324 seconds (Total)
#> Chain 2: 
#> Warning: There were 6 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit2)
#> Warning: There were 6 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ gp(x0) + x1 + gp(x2, k = 10) + x3 
#>    Data: dat (Number of observations: 30) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Gaussian Process Hyperparameters:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sdgp(gpx0)       0.82      0.65     0.03     2.41 1.00      885     1113
#> sdgp(gpx2)       5.53      2.24     2.43    10.85 1.00     1226     1404
#> lscale(gpx0)     0.23      3.53     0.01     0.61 1.00      991      561
#> lscale(gpx2)     0.07      0.05     0.01     0.20 1.00     1749     1221
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     6.27      2.21     1.65    10.39 1.00     1566     1177
#> x1            3.25      1.71    -0.08     6.75 1.00     2612     1391
#> x3           -1.08      1.66    -4.49     1.99 1.00     2000     1245
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     1.93      0.35     1.36     2.70 1.00     1045     1222
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
me2 <- conditional_effects(fit2, ndraws = 200, spaghetti = TRUE)
plot(me2, ask = FALSE, points = TRUE)





# fit a multivariate GP model with Matern 3/2 kernel
fit3 <- brm(y ~ gp(x1, x2, cov = "matern32"), dat, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000102 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.02 seconds.
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
#> Chain 1:  Elapsed Time: 3.675 seconds (Warm-up)
#> Chain 1:                2.471 seconds (Sampling)
#> Chain 1:                6.146 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 8.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.84 seconds.
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
#> Chain 2:  Elapsed Time: 3.69 seconds (Warm-up)
#> Chain 2:                3.928 seconds (Sampling)
#> Chain 2:                7.618 seconds (Total)
#> Chain 2: 
#> Warning: There were 61 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
summary(fit3)
#> Warning: There were 61 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ gp(x1, x2, cov = "matern32") 
#>    Data: dat (Number of observations: 30) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Gaussian Process Hyperparameters:
#>                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sdgp(gpx1x2)       3.82      1.16     2.15     6.74 1.01      464      436
#> lscale(gpx1x2)     0.28      0.18     0.12     0.69 1.01      310      533
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     7.57      1.82     3.76    11.14 1.00      543     1029
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     1.62      0.47     0.80     2.58 1.02      104       93
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
me3 <- conditional_effects(fit3, ndraws = 200, spaghetti = TRUE)
plot(me3, ask = FALSE, points = TRUE)




# compare model fit
loo(fit1, fit2, fit3)
#> Warning: Found 4 observations with a pareto_k > 0.7 in model 'fit1'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> Warning: Found 2 observations with a pareto_k > 0.7 in model 'fit2'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> Warning: Found 8 observations with a pareto_k > 0.7 in model 'fit3'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations. 
#> Output of model 'fit1':
#> 
#> Computed from 2000 by 30 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -70.8 4.8
#> p_loo        10.5 3.6
#> looic       141.6 9.6
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.0]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     26    86.7%   187     
#>    (0.7, 1]   (bad)       3    10.0%   <NA>    
#>    (1, Inf)   (very bad)  1     3.3%   <NA>    
#> See help('pareto-k-diagnostic') for details.
#> 
#> Output of model 'fit2':
#> 
#> Computed from 2000 by 30 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -67.2 3.0
#> p_loo         9.2 1.5
#> looic       134.4 6.0
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 1.1]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     28    93.3%   168     
#>    (0.7, 1]   (bad)       2     6.7%   <NA>    
#>    (1, Inf)   (very bad)  0     0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.
#> 
#> Output of model 'fit3':
#> 
#> Computed from 2000 by 30 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -65.9 3.4
#> p_loo        15.8 2.6
#> looic       131.8 6.8
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 0.4]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     22    73.3%   10      
#>    (0.7, 1]   (bad)       8    26.7%   <NA>    
#>    (1, Inf)   (very bad)  0     0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.
#> 
#> Model comparisons:
#>  model elpd_diff se_diff p_worse diag_diff      diag_elpd
#>   fit3       0.0     0.0      NA           8 k_psis > 0.7
#>   fit2      -1.3     2.6    0.69   N < 100 2 k_psis > 0.7
#>   fit1      -4.9     3.1    0.94   N < 100 4 k_psis > 0.7
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# simulate data with a factor covariate
dat2 <- mgcv::gamSim(4, n = 90, scale = 2)
#> Factor `by' variable example

# fit separate gaussian processes for different levels of 'fac'
fit4 <- brm(y ~ gp(x2, by = fac), dat2, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000146 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.46 seconds.
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
#> Chain 1:  Elapsed Time: 2.838 seconds (Warm-up)
#> Chain 1:                2.672 seconds (Sampling)
#> Chain 1:                5.51 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.0001 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1 seconds.
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
#> Chain 2:  Elapsed Time: 2.912 seconds (Warm-up)
#> Chain 2:                2.687 seconds (Sampling)
#> Chain 2:                5.599 seconds (Total)
#> Chain 2: 
#> Warning: There were 7 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
summary(fit4)
#> Warning: There were 7 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: y ~ gp(x2, by = fac) 
#>    Data: dat2 (Number of observations: 90) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Gaussian Process Hyperparameters:
#>                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sdgp(gpx2fac1)       0.83      0.56     0.03     2.04 1.00      529      787
#> sdgp(gpx2fac2)       2.43      0.96     0.98     4.59 1.00      979     1010
#> sdgp(gpx2fac3)       2.86      1.12     1.32     5.42 1.00      867      933
#> lscale(gpx2fac1)     0.07      0.17     0.01     0.41 1.00      458      642
#> lscale(gpx2fac2)     0.11      0.47     0.02     0.40 1.00      471      507
#> lscale(gpx2fac3)     0.14      0.08     0.05     0.37 1.00      496      460
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept     1.25      0.42     0.41     2.08 1.00     1409     1209
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     2.00      0.21     1.59     2.40 1.00      628      773
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_effects(fit4), points = TRUE)



# }
```
