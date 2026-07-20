# GPU support in Stan via OpenCL

Use OpenCL for GPU support in Stan via the brms interface. Only some
Stan functions can be run on a GPU at this point and so a lot of brms
models won't benefit from OpenCL for now.

## Usage

``` r
opencl(ids = NULL)
```

## Arguments

- ids:

  (integer vector of length 2) The platform and device IDs of the OpenCL
  device to use for fitting. If you don't know the IDs of your OpenCL
  device, `c(0,0)` is most likely what you need.

## Value

A `brmsopencl` object which can be passed to the `opencl` argument of
`brm` and related functions.

## Details

For more details on OpenCL in Stan, check out
<https://mc-stan.org/docs/2_26/cmdstan-guide/parallelization.html#opencl>
as well as <https://mc-stan.org/docs/2_26/stan-users-guide/opencl.html>.

## Examples

``` r
# \dontrun{
# this model just serves as an illustration
# OpenCL may not actually speed things up here
fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
           data = epilepsy, family = poisson(),
           chains = 2, cores = 2, opencl = opencl(c(0, 0)),
           backend = "cmdstanr")
#> /usr/bin/ld: cannot find -lOpenCL: No such file or directory
#> collect2: error: ld returned 1 exit status
#> make: *** [make/program:82: /tmp/RtmpBKJ0Ra/model-25ba98b3769] Error 1
#> Error: An error occured during compilation! See the message above for more information.
summary(fit)
#> Error: object 'fit' not found
# }
```
