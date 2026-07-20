# Posterior Probabilities of Mixture Component Memberships

Compute the posterior probabilities of mixture component memberships for
each observation including uncertainty estimates.

## Usage

``` r
# S3 method for class 'brmsfit'
pp_mixture(
  x,
  newdata = NULL,
  re_formula = NULL,
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  log = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)

pp_mixture(x, ...)
```

## Arguments

- x:

  An R object usually of class `brmsfit`.

- newdata:

  An optional data.frame for which to evaluate predictions. If `NULL`
  (default), the original data of the model is used. `NA` values within
  factors (excluding grouping variables) are interpreted as if all dummy
  variables of this factor are zero. This allows, for instance, to make
  predictions of the grand mean when using sum coding. `NA` values
  within grouping variables are treated as a new level.

- re_formula:

  formula containing group-level effects to be considered in the
  prediction. If `NULL` (default), include all group-level effects; if
  `NA` or `~0`, include no group-level effects.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- ndraws:

  Positive integer indicating how many posterior draws should be used.
  If `NULL` (the default) all draws are used. Ignored if `draw_ids` is
  not `NULL`.

- draw_ids:

  An integer vector specifying the posterior draws to be used. If `NULL`
  (the default), all draws are used.

- log:

  Logical; Indicates whether to return probabilities on the log-scale.

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead. Only used if `summary` is `TRUE`.

- probs:

  The percentiles to be computed by the `quantile` function. Only used
  if `summary` is `TRUE`.

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

If `summary = TRUE`, an N x E x K array, where N is the number of
observations, K is the number of mixture components, and E is equal to
`length(probs) + 2`. If `summary = FALSE`, an S x N x K array, where S
is the number of posterior draws.

## Details

The returned probabilities can be written as \\P(K_n = k \| Y_n)\\, that
is the posterior probability that observation n originates from
component k. They are computed using Bayes' Theorem \$\$P(K_n = k \|
Y_n) = P(Y_n \| K_n = k) P(K_n = k) / P(Y_n),\$\$ where \\P(Y_n \| K_n =
k)\\ is the (posterior) likelihood of observation n for component k,
\\P(K_n = k)\\ is the (posterior) mixing probability of component k
(i.e. parameter `theta<k>`), and \$\$P(Y_n) = \sum\_{k=1,...,K} P(Y_n \|
K_n = k) P(K_n = k)\$\$ is a normalizing constant.

## Examples

``` r
# \dontrun{
## simulate some data
set.seed(1234)
dat <- data.frame(
  y = c(rnorm(100), rnorm(50, 2)),
  x = rnorm(150)
)
## fit a simple normal mixture model
mix <- mixture(gaussian, nmix = 2)
#> Setting order = 'mu' for mixtures of the same family.
prior <- c(
  prior(normal(0, 5), Intercept, nlpar = mu1),
  prior(normal(0, 5), Intercept, nlpar = mu2),
  prior(dirichlet(2, 2), theta)
)
fit1 <- brm(bf(y ~ x), dat, family = mix,
            prior = prior, chains = 2, init = 0)
#> Error in .validate_prior(prior, bframe = bframe, sample_prior = sample_prior): The following priors do not correspond to any model parameter: 
#> Intercept_mu1 ~ normal(0, 5)
#> Intercept_mu2 ~ normal(0, 5)
#> Function 'default_prior' might be helpful to you.
summary(fit1)
#> Error: object 'fit1' not found

## compute the membership probabilities
ppm <- pp_mixture(fit1)
#> Error: object 'fit1' not found
str(ppm)
#> Error: object 'ppm' not found

## extract point estimates for each observation
head(ppm[, 1, ])
#> Error: object 'ppm' not found

## classify every observation according to
## the most likely component
apply(ppm[, 1, ], 1, which.max)
#> Error: object 'ppm' not found
# }
```
