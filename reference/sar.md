# Spatial simultaneous autoregressive (SAR) structures

Set up an spatial simultaneous autoregressive (SAR) term in brms. The
function does not evaluate its arguments – it exists purely to help set
up a model with SAR terms.

## Usage

``` r
sar(M, type = "lag")
```

## Arguments

- M:

  An object specifying the spatial weighting matrix. Can be either the
  spatial weight matrix itself or an object of class `listw` or `nb`,
  from which the spatial weighting matrix can be computed.

- type:

  Type of the SAR structure. Either `"lag"` (for SAR of the response
  values) or `"error"` (for SAR of the residuals). More information is
  provided in the 'Details' section.

## Value

An object of class `'sar_term'`, which is a list of arguments to be
interpreted by the formula parsing functions of brms.

## Details

The `lagsar` structure implements SAR of the response values: \$\$y =
\rho W y + \eta + e\$\$ The `errorsar` structure implements SAR of the
residuals: \$\$y = \eta + u, u = \rho W u + e\$\$ In the above
equations, \\\eta\\ is the predictor term and \\e\\ are independent
normally or t-distributed residuals. Currently, only families `gaussian`
and `student` support SAR structures.

## See also

[`autocor-terms`](https://paulbuerkner.com/brms/reference/autocor-terms.md)

## Examples

``` r
# \dontrun{
data(oldcol, package = "spdep")
fit1 <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "lag"),
            data = COL.OLD, data2 = list(COL.nb = COL.nb),
            chains = 2, cores = 2)
#> Compiling Stan program...
#> Start sampling
summary(fit1)
plot(fit1)


fit2 <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "error"),
            data = COL.OLD, data2 = list(COL.nb = COL.nb),
            chains = 2, cores = 2)
#> Compiling Stan program...
#> Start sampling
summary(fit2)
#> 1333     1201
#> HOVAL        -0.30      0.10    -0.48    -0.11 1.00     1753     1481
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma    10.39      1.13     8.52    12.93 1.00     2034     1341
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
#> ng)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.473 seconds (Warm-up)
#> Chain 1:                0.294 seconds (Sampling)
#> Chain 1:                0.767 seconds (Total)
#> Chain 1: 
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.474 seconds (Warm-up)
#> Chain 2:                0.306 seconds (Sampling)
#> Chain 2:                0.78 seconds (Total)
#> Chain 2: 
plot(fit2)

# }
```
