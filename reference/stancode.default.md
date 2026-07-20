# Stan Code for brms Models

Generate Stan code for brms models

## Usage

``` r
# Default S3 method
stancode(
  object,
  data,
  family = gaussian(),
  prior = NULL,
  autocor = NULL,
  data2 = NULL,
  cov_ranef = NULL,
  sparse = NULL,
  sample_prior = "no",
  stanvars = NULL,
  stan_funs = NULL,
  knots = NULL,
  drop_unused_levels = TRUE,
  threads = getOption("brms.threads", NULL),
  normalize = getOption("brms.normalize", TRUE),
  save_model = NULL,
  ...
)
```

## Arguments

- object:

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

- prior:

  One or more `brmsprior` objects created by
  [`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md) or
  related functions and combined using the `c` method or the `+`
  operator. See also
  [`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md)
  for more help.

- autocor:

  (Deprecated) An optional
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md)
  object describing the correlation structure within the response
  variable (i.e., the 'autocorrelation'). See the documentation of
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md) for
  a description of the available correlation structures. Defaults to
  `NULL`, corresponding to no correlations. In multivariate models,
  `autocor` might also be a list of autocorrelation structures. It is
  now recommend to specify autocorrelation terms directly within
  `formula`. See
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  for more details.

- data2:

  A named `list` of objects containing data, which cannot be passed via
  argument `data`. Required for some objects used in autocorrelation
  structures to specify dependency structures as well as for
  within-group covariance matrices.

- cov_ranef:

  (Deprecated) A list of matrices that are proportional to the (within)
  covariance structure of the group-level effects. The names of the
  matrices should correspond to columns in `data` that are used as
  grouping factors. All levels of the grouping factor should appear as
  rownames of the corresponding matrix. This argument can be used, among
  others to model pedigrees and phylogenetic effects. It is now
  recommended to specify those matrices in the formula interface using
  the [`gr`](https://paulbuerkner.com/brms/reference/gr.md) and related
  functions. See
  [`vignette("brms_phylogenetics")`](https://paulbuerkner.com/brms/articles/brms_phylogenetics.md)
  for more details.

- sparse:

  (Deprecated) Logical; indicates whether the population-level design
  matrices should be treated as sparse (defaults to `FALSE`). For design
  matrices with many zeros, this can considerably reduce required
  memory. Sampling speed is currently not improved or even slightly
  decreased. It is now recommended to use the `sparse` argument of
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  and related functions.

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

- stanvars:

  An optional `stanvars` object generated by function
  [`stanvar`](https://paulbuerkner.com/brms/reference/stanvar.md) to
  define additional variables for use in Stan's program blocks.

- stan_funs:

  (Deprecated) An optional character string containing self-defined Stan
  functions, which will be included in the functions block of the
  generated Stan code. It is now recommended to use the `stanvars`
  argument for this purpose instead.

- knots:

  Optional list containing user specified knot values to be used for
  basis construction of smoothing terms. See
  [`gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html) for more details.

- drop_unused_levels:

  Should unused factors levels in the data be dropped? Defaults to
  `TRUE`.

- threads:

  Number of threads to use in within-chain parallelization. For more
  control over the threading process, `threads` may also be a
  `brmsthreads` object created by
  [`threading`](https://paulbuerkner.com/brms/reference/threading.md).
  Within-chain parallelization is experimental! We recommend its use
  only if you are experienced with Stan's `reduce_sum` function and have
  a slow running model that cannot be sped up by any other means. Can be
  set globally for the current R session via the `"brms.threads"` option
  (see [`options`](https://rdrr.io/r/base/options.html)).

- normalize:

  Logical. Indicates whether normalization constants should be included
  in the Stan code (defaults to `TRUE`). Setting it to `FALSE` requires
  Stan version \>= 2.25 to work. If `FALSE`, sampling efficiency may be
  increased but some post processing functions such as
  [`bridge_sampler`](https://paulbuerkner.com/brms/reference/bridge_sampler.brmsfit.md)
  will not be available. Can be controlled globally for the current R
  session via the \`brms.normalize\` option.

- save_model:

  Either `NULL` or a character string. In the latter case, the model's
  Stan code is saved via
  [`cat`](https://paulbuerkner.com/brms/reference/addition-terms.md) in
  a text file named after the string supplied in `save_model`.

- ...:

  Other arguments for internal usage only.

## Value

A character string containing the fully commented Stan code to fit a
brms model. It is of class `c("character", "brmsmodel")` to facilitate
pretty printing.

## Examples

``` r
stancode(rating ~ treat + period + carry + (1|subject),
         data = inhaler, family = "cumulative")
#> // generated with brms 2.23.1
#> functions {
#>   /* cumulative-logit log-PDF for a single response
#>    * Args:
#>    *   y: response category
#>    *   mu: latent mean parameter
#>    *   disc: discrimination parameter
#>    *   thres: ordinal thresholds
#>    * Returns:
#>    *   a scalar to be added to the log posterior
#>    */
#>    real cumulative_logit_lpmf(int y, real mu, real disc, vector thres) {
#>      int nthres = num_elements(thres);
#>      if (y == 1) {
#>        return log_inv_logit(disc * (thres[1] - mu));
#>      } else if (y == nthres + 1) {
#>        return log1m_inv_logit(disc * (thres[nthres] - mu));
#>      } else {
#>        return log_inv_logit_diff(disc * (thres[y] - mu), disc * (thres[y - 1] - mu));
#>      }
#>    }
#>   /* cumulative-logit log-PDF for a single response and merged thresholds
#>    * Args:
#>    *   y: response category
#>    *   mu: latent mean parameter
#>    *   disc: discrimination parameter
#>    *   thres: vector of merged ordinal thresholds
#>    *   j: start and end index for the applid threshold within 'thres'
#>    * Returns:
#>    *   a scalar to be added to the log posterior
#>    */
#>    real cumulative_logit_merged_lpmf(int y, real mu, real disc, vector thres, array[] int j) {
#>      return cumulative_logit_lpmf(y | mu, disc, thres[j[1]:j[2]]);
#>    }
#>   /* ordered-logistic log-PDF for a single response and merged thresholds
#>    * Args:
#>    *   y: response category
#>    *   mu: latent mean parameter
#>    *   thres: vector of merged ordinal thresholds
#>    *   j: start and end index for the applid threshold within 'thres'
#>    * Returns:
#>    *   a scalar to be added to the log posterior
#>    */
#>    real ordered_logistic_merged_lpmf(int y, real mu, vector thres, array[] int j) {
#>      return ordered_logistic_lpmf(y | mu, thres[j[1]:j[2]]);
#>    }
#> 
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   array[N] int Y;  // response variable
#>   int<lower=2> nthres;  // number of thresholds
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   // data for group-level effects of ID 1
#>   int<lower=1> N_1;  // number of grouping levels
#>   int<lower=1> M_1;  // number of coefficients per level
#>   array[N] int<lower=1> J_1;  // grouping indicator per observation
#>   // group-level predictor values
#>   vector[N] Z_1_1;
#>   int prior_only;  // should the likelihood be ignored?
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 1:K) {
#>     means_X[i] = mean(X[, i]);
#>     Xc[, i] = X[, i] - means_X[i];
#>   }
#> }
#> parameters {
#>   vector[Kc] b;  // regression coefficients
#>   ordered[nthres] Intercept;  // temporary thresholds for centered predictors
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#> }
#> transformed parameters {
#>   real disc = 1;  // discrimination parameters
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
#>   lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
#>     - 1 * student_t_lccdf(0 | 3, 0, 2.5);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Xc * b;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#>     }
#>     for (n in 1:N) {
#>       target += ordered_logistic_lpmf(Y[n] | mu[n], Intercept);
#>     }
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#> }
#> generated quantities {
#>   // compute actual thresholds
#>   vector[nthres] b_Intercept = Intercept + dot_product(means_X, b);
#> }

stancode(count ~ zAge + zBase * Trt + (1|patient),
         data = epilepsy, family = "poisson")
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   array[N] int Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   // data for group-level effects of ID 1
#>   int<lower=1> N_1;  // number of grouping levels
#>   int<lower=1> M_1;  // number of coefficients per level
#>   array[N] int<lower=1> J_1;  // grouping indicator per observation
#>   // group-level predictor values
#>   vector[N] Z_1_1;
#>   int prior_only;  // should the likelihood be ignored?
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X without an intercept
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 2:K) {
#>     means_X[i - 1] = mean(X[, i]);
#>     Xc[, i - 1] = X[, i] - means_X[i - 1];
#>   }
#> }
#> parameters {
#>   vector[Kc] b;  // regression coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#> }
#> transformed parameters {
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   lprior += student_t_lpdf(Intercept | 3, 1.4, 2.5);
#>   lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
#>     - 1 * student_t_lccdf(0 | 3, 0, 2.5);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Intercept;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#>     }
#>     target += poisson_log_glm_lpmf(Y | Xc, mu, b);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }
```
