# Default priors for Bayesian models

`default_prior` is a generic function that can be used to get default
priors for Bayesian models. Its original use is within the brms package,
but new methods for use with objects from other packages can be
registered to the same generic.

## Usage

``` r
default_prior(object, ...)

get_prior(formula, ...)
```

## Arguments

- object:

  An object whose class will determine which method will be used. A
  symbolic description of the model to be fitted.

- ...:

  Further arguments passed to the specific method.

- formula:

  Synonym of `object` for use in `get_prior`.

## Value

Usually, a `brmsprior` object. See
[`default_prior.default`](https://paulbuerkner.com/brms/reference/default_prior.default.md)
for more details.

## Details

See
[`default_prior.default`](https://paulbuerkner.com/brms/reference/default_prior.default.md)
for the default method applied for brms models. You can view the
available methods by typing `methods(default_prior)`.

## See also

[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md),
[`default_prior.default`](https://paulbuerkner.com/brms/reference/default_prior.default.md)

## Examples

``` r
## get all parameters and parameters classes to define priors on
(prior <- default_prior(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
                        data = epilepsy, family = poisson()))
#>                   prior     class       coef   group resp dpar nlpar lb ub tag
#>  student_t(3, 1.4, 2.5) Intercept                                             
#>                  (flat)         b                                             
#>                  (flat)         b       Trt1                                  
#>                  (flat)         b       zAge                                  
#>                  (flat)         b      zBase                                  
#>                  (flat)         b zBase:Trt1                                  
#>    student_t(3, 0, 2.5)        sd                                     0       
#>    student_t(3, 0, 2.5)        sd                obs                  0       
#>    student_t(3, 0, 2.5)        sd  Intercept     obs                  0       
#>    student_t(3, 0, 2.5)        sd            patient                  0       
#>    student_t(3, 0, 2.5)        sd  Intercept patient                  0       
#>        source
#>       default
#>       default
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>       default
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
```
