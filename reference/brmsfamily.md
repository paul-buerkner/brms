# Special Family Functions for brms Models

Family objects provide a convenient way to specify the details of the
models used by many model fitting functions. The family functions
presented here are for use with brms only and will \*\*not\*\* work with
other model fitting functions such as `glm` or `glmer`. However, the
standard family functions as described in
[`family`](https://rdrr.io/r/stats/family.html) will work with brms. You
can also specify custom families for use in brms with the
[`custom_family`](https://paulbuerkner.com/brms/reference/custom_family.md)
function.

## Usage

``` r
brmsfamily(
  family,
  link = NULL,
  link_sigma = "log",
  link_shape = "log",
  link_nu = "logm1",
  link_phi = "log",
  link_kappa = "log",
  link_beta = "log",
  link_zi = "logit",
  link_hu = "logit",
  link_zoi = "logit",
  link_coi = "logit",
  link_disc = "log",
  link_bs = "log",
  link_ndt = "log",
  link_bias = "logit",
  link_xi = "log1p",
  link_alpha = "identity",
  link_quantile = "logit",
  threshold = "flexible",
  refcat = NULL
)

student(link = "identity", link_sigma = "log", link_nu = "logm1")

bernoulli(link = "logit")

beta_binomial(link = "logit", link_phi = "log")

negbinomial(link = "log", link_shape = "log")

geometric(link = "log")

lognormal(link = "identity", link_sigma = "log")

shifted_lognormal(link = "identity", link_sigma = "log", link_ndt = "log")

skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity")

exponential(link = "log")

weibull(link = "log", link_shape = "log")

frechet(link = "log", link_nu = "logm1")

gen_extreme_value(link = "identity", link_sigma = "log", link_xi = "log1p")

exgaussian(link = "identity", link_sigma = "log", link_beta = "log")

wiener(
  link = "identity",
  link_bs = "log",
  link_ndt = "log",
  link_bias = "logit"
)

Beta(link = "logit", link_phi = "log")

xbeta(link = "logit", link_phi = "log", link_kappa = "log")

dirichlet(link = "logit", link_phi = "log", refcat = NULL)

logistic_normal(link = "identity", link_sigma = "log", refcat = NULL)

von_mises(link = "tan_half", link_kappa = "log")

asym_laplace(link = "identity", link_sigma = "log", link_quantile = "logit")

cox(link = "log")

hurdle_poisson(link = "log", link_hu = "logit")

hurdle_negbinomial(link = "log", link_shape = "log", link_hu = "logit")

hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit")

hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit")

hurdle_cumulative(
  link = "logit",
  link_hu = "logit",
  link_disc = "log",
  threshold = "flexible"
)

zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit")

zero_one_inflated_beta(
  link = "logit",
  link_phi = "log",
  link_zoi = "logit",
  link_coi = "logit"
)

zero_inflated_poisson(link = "log", link_zi = "logit")

zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

zero_inflated_binomial(link = "logit", link_zi = "logit")

zero_inflated_beta_binomial(
  link = "logit",
  link_phi = "log",
  link_zi = "logit"
)

categorical(link = "logit", refcat = NULL)

multinomial(link = "logit", refcat = NULL)

dirichlet_multinomial(link = "logit", link_phi = "log", refcat = NULL)

cumulative(link = "logit", link_disc = "log", threshold = "flexible")

sratio(link = "logit", link_disc = "log", threshold = "flexible")

cratio(link = "logit", link_disc = "log", threshold = "flexible")

acat(link = "logit", link_disc = "log", threshold = "flexible")
```

## Arguments

- family:

  A character string naming the distribution family of the response
  variable to be used in the model. Currently, the following families
  are supported: `gaussian`, `student`, `binomial`, `bernoulli`,
  `beta-binomial`, `poisson`, `negbinomial`, `geometric`, `Gamma`,
  `skew_normal`, `lognormal`, `shifted_lognormal`, `exgaussian`,
  `wiener`, `inverse.gaussian`, `exponential`, `weibull`, `frechet`,
  `Beta`, `dirichlet`, `von_mises`, `asym_laplace`, `gen_extreme_value`,
  `categorical`, `multinomial`, `dirichlet_multinomial`, `cumulative`,
  `cratio`, `sratio`, `acat`, `hurdle_poisson`, `hurdle_negbinomial`,
  `hurdle_gamma`, `hurdle_lognormal`, `hurdle_cumulative`,
  `zero_inflated_binomial`, `zero_inflated_beta_binomial`,
  `zero_inflated_beta`, `zero_inflated_negbinomial`,
  `zero_inflated_poisson`, `zero_one_inflated_beta`, and `xbeta`.

- link:

  A specification for the model link function. This can be a
  name/expression or character string. See the 'Details' section for
  more information on link functions supported by each family.

- link_sigma:

  Link of auxiliary parameter `sigma` if being predicted.

- link_shape:

  Link of auxiliary parameter `shape` if being predicted.

- link_nu:

  Link of auxiliary parameter `nu` if being predicted.

- link_phi:

  Link of auxiliary parameter `phi` if being predicted.

- link_kappa:

  Link of auxiliary parameter `kappa` if being predicted.

- link_beta:

  Link of auxiliary parameter `beta` if being predicted.

- link_zi:

  Link of auxiliary parameter `zi` if being predicted.

- link_hu:

  Link of auxiliary parameter `hu` if being predicted.

- link_zoi:

  Link of auxiliary parameter `zoi` if being predicted.

- link_coi:

  Link of auxiliary parameter `coi` if being predicted.

- link_disc:

  Link of auxiliary parameter `disc` if being predicted.

- link_bs:

  Link of auxiliary parameter `bs` if being predicted.

- link_ndt:

  Link of auxiliary parameter `ndt` if being predicted.

- link_bias:

  Link of auxiliary parameter `bias` if being predicted.

- link_xi:

  Link of auxiliary parameter `xi` if being predicted.

- link_alpha:

  Link of auxiliary parameter `alpha` if being predicted.

- link_quantile:

  Link of auxiliary parameter `quantile` if being predicted.

- threshold:

  A character string indicating the type of thresholds (i.e. intercepts)
  used in an ordinal model. `"flexible"` provides the standard
  unstructured thresholds, `"equidistant"` restricts the distance
  between consecutive thresholds to the same value, and `"sum_to_zero"`
  ensures the thresholds sum to zero.

- refcat:

  Optional name of the reference response category used in
  `categorical`, `multinomial`, `dirichlet`, `dirichlet_multinomial` and
  `logistic_normal` models. If `NULL` (the default), the first category
  is used as the reference. If `NA`, all categories will be predicted,
  which requires strong priors or carefully specified predictor terms in
  order to lead to an identified model.

## Details

Below, we list common use cases for the different families. This list is
not ment to be exhaustive.

- Family `gaussian` can be used for linear regression.

- Family `student` can be used for robust linear regression that is less
  influenced by outliers.

- Family `skew_normal` can handle skewed responses in linear regression.

- Families `poisson`, `negbinomial`, and `geometric` can be used for
  regression of unbounded count data.

- Families `bernoulli`, `binomial`, and `beta_binomial` can be used for
  binary regression (i.e., most commonly logistic regression).

- Families `categorical`, `multinomial` and `dirichlet_multinomial` can
  be used for multi-logistic regression when there are more than two
  possible outcomes.

- Families `cumulative`, `cratio` ('continuation ratio'), `sratio`
  ('stopping ratio'), and `acat` ('adjacent category') leads to ordinal
  regression.

- Families `Gamma`, `weibull`, `exponential`, `lognormal`, `frechet`,
  `inverse.gaussian`, and `cox` (Cox proportional hazards model) can be
  used (among others) for time-to-event regression also known as
  survival regression.

- Families `weibull`, `frechet`, and `gen_extreme_value` ('generalized
  extreme value') allow for modeling extremes.

- Families `beta`, `dirichlet`, and `logistic_normal` can be used to
  model responses representing rates or probabilities.

- Family `xbeta` extends the `beta` family to support `[0, 1]` responses
  with exact `0`s and / or `1`s, when each response takes values `0`,
  `1`, and `(0, 1)` according to a single process. If there is merit in
  assuming that 0 and 1 values arise from different processes than
  `(0, 1)` values, then the `zero_inflated_beta`,
  `zero_one_inflated_beta` families provide more flexibility. For
  details see Kosmidis & Zeileis (2024).

- Family `asym_laplace` allows for quantile regression when fixing the
  auxiliary `quantile` parameter to the quantile of interest.

- Family `exgaussian` ('exponentially modified Gaussian') and
  `shifted_lognormal` are especially suited to model reaction times.

- Family `wiener` provides an implementation of the Wiener diffusion
  model. For this family, the main formula predicts the drift parameter
  'delta' and all other parameters are modeled as auxiliary parameters
  (see
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  for details).

- Families `hurdle_poisson`, `hurdle_negbinomial`, `hurdle_gamma`,
  `hurdle_lognormal`, `zero_inflated_poisson`,
  `zero_inflated_negbinomial`, `zero_inflated_binomial`,
  `zero_inflated_beta_binomial`, `zero_inflated_beta`,
  `zero_one_inflated_beta`, and `hurdle_cumulative` allow to estimate
  zero-inflated and hurdle models. These models can be very helpful when
  there are many zeros in the data (or ones in case of one-inflated
  models) that cannot be explained by the primary distribution of the
  response.

Below, we list all possible links for each family. The first link
mentioned for each family is the default.

- Families `gaussian`, `student`, `skew_normal`, `exgaussian`,
  `asym_laplace`, and `gen_extreme_value` support the links (as names)
  `identity`, `log`, `inverse`, and `softplus`.

- Families `poisson`, `negbinomial`, `geometric`,
  `zero_inflated_poisson`, `zero_inflated_negbinomial`,
  `hurdle_poisson`, and `hurdle_negbinomial` support `log`, `identity`,
  `sqrt`, and `softplus`.

- Families `binomial`, `bernoulli`, `beta_binomial`,
  `zero_inflated_binomial`, `zero_inflated_beta_binomial`, `Beta`,
  `zero_inflated_beta`, `zero_one_inflated_beta`, and `xbeta` support
  `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `identity`,
  and `log`.

- Families `cumulative`, `cratio`, `sratio`, `acat`, and
  `hurdle_cumulative` support `logit`, `probit`, `probit_approx`,
  `cloglog`, and `cauchit`.

- Families `categorical`, `multinomial`, `dirichlet_multinomial` and
  `dirichlet` support `logit`.

- Families `Gamma`, `weibull`, `exponential`, `frechet`, and
  `hurdle_gamma` support `log`, `identity`, `inverse`, and `softplus`.

- Families `lognormal` and `hurdle_lognormal` support `identity` and
  `inverse`.

- Family `logistic_normal` supports `identity`.

- Family `inverse.gaussian` supports `1/mu^2`, `inverse`, `identity`,
  `log`, and `softplus`.

- Family `von_mises` supports `tan_half` and `identity`.

- Family `cox` supports `log`, `identity`, and `softplus` for the
  proportional hazards parameter.

- Family `wiener` supports `identity`, `log`, and `softplus` for the
  main parameter which represents the drift rate.

Please note that when calling the
[`Gamma`](https://rdrr.io/r/stats/family.html) family function of the
stats package, the default link will be `inverse` instead of `log`
although the latter is the default in brms. Also, when using the family
functions `gaussian`, `binomial`, `poisson`, and `Gamma` of the stats
package (see [`family`](https://rdrr.io/r/stats/family.html)), special
link functions such as `softplus` or `cauchit` won't work. In this case,
you have to use `brmsfamily` to specify the family with corresponding
link function.

## References

Kosmidis I, Zeileis A (2024). Extended-Support Beta Regression for \[0,
1\] Responses. *arXiv Preprint*.
[doi:10.48550/arXiv.2409.07233](https://doi.org/10.48550/arXiv.2409.07233)

## See also

[`brm`](https://paulbuerkner.com/brms/reference/brm.md),
[`family`](https://rdrr.io/r/stats/family.html),
[`customfamily`](https://paulbuerkner.com/brms/reference/custom_family.md)

## Examples

``` r
 # create a family object
 (fam1 <- student("log"))
#> 
#> Family: student 
#> Link function: log 
#> 
 # alternatively use the brmsfamily function
 (fam2 <- brmsfamily("student", "log"))
#> 
#> Family: student 
#> Link function: log 
#> 
 # both leads to the same object
 identical(fam1, fam2)
#> [1] FALSE
```
