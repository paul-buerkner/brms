# Validate Prior for brms Models

Validate priors supplied by the user. Return a complete set of priors
for the given model, including default priors.

## Usage

``` r
validate_prior(
  prior,
  formula,
  data,
  family = gaussian(),
  sample_prior = "no",
  data2 = NULL,
  knots = NULL,
  drop_unused_levels = TRUE,
  ...
)
```

## Arguments

- prior:

  One or more `brmsprior` objects created by
  [`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md) or
  related functions and combined using the `c` method or the `+`
  operator. See also
  [`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md)
  for more help.

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html),
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
  or
  [`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)
  (or one that can be coerced to that classes): A symbolic description
  of the model to be fitted. The details of model specification are
  explained in
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

- data:

  An object of class `data.frame` (or one that can be coerced to that
  class) containing data of all variables used in the model.

- family:

  A description of the response distribution and link function to be
  used in the model. This can be a family function, a call to a family
  function or a character string naming the family. Every family
  function has a `link` argument allowing to specify the link function
  to be applied on the response variable. If not specified, default
  links are used. For details of supported families see
  [`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).
  By default, a linear `gaussian` model is applied. In multivariate
  models, `family` might also be a list of families.

- sample_prior:

  Indicate if draws from priors should be drawn additionally to the
  posterior draws. Options are `"no"` (the default), `"yes"`, and
  `"only"`. Among others, these draws can be used to calculate Bayes
  factors for point hypotheses via
  [`hypothesis`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.md).
  Please note that improper priors are not sampled, including the
  default improper priors used by `brm`. See
  [`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md) on
  how to set (proper) priors. Please also note that prior draws for the
  overall intercept are not obtained by default for technical reasons.
  See
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  how to obtain prior draws for the intercept. If `sample_prior` is set
  to `"only"`, draws are drawn solely from the priors ignoring the
  likelihood, which allows among others to generate draws from the prior
  predictive distribution. In this case, all parameters must have proper
  priors.

- data2:

  A named `list` of objects containing data, which cannot be passed via
  argument `data`. Required for some objects used in autocorrelation
  structures to specify dependency structures as well as for
  within-group covariance matrices.

- knots:

  Optional list containing user specified knot values to be used for
  basis construction of smoothing terms. See
  [`gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html) for more details.

- drop_unused_levels:

  Should unused factors levels in the data be dropped? Defaults to
  `TRUE`.

- ...:

  Other arguments for internal usage only.

## Value

An object of class `brmsprior`.

## See also

[`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md),
[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md).

## Examples

``` r
prior1 <- prior(normal(0,10), class = b) +
  prior(cauchy(0,2), class = sd)
validate_prior(prior1, count ~ zAge + zBase * Trt + (1|patient),
               data = epilepsy, family = poisson())
#>                   prior     class       coef   group resp dpar nlpar lb ub tag
#>  student_t(3, 1.4, 2.5) Intercept                                             
#>           normal(0, 10)         b                                             
#>           normal(0, 10)         b       Trt1                                  
#>           normal(0, 10)         b       zAge                                  
#>           normal(0, 10)         b      zBase                                  
#>           normal(0, 10)         b zBase:Trt1                                  
#>            cauchy(0, 2)        sd                                     0       
#>            cauchy(0, 2)        sd            patient                  0       
#>            cauchy(0, 2)        sd  Intercept patient                  0       
#>        source
#>       default
#>          user
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>          user
#>  (vectorized)
#>  (vectorized)
```
