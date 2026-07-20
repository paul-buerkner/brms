# Set up a model formula for use in brms

Set up a model formula for use in the brms package allowing to define
(potentially non-linear) additive multilevel models for all parameters
of the assumed response distribution.

## Usage

``` r
brmsformula(
  formula,
  ...,
  flist = NULL,
  family = NULL,
  autocor = NULL,
  nl = NULL,
  loop = NULL,
  center = NULL,
  cmc = NULL,
  sparse = NULL,
  decomp = NULL,
  unused = NULL
)
```

## Arguments

- formula:

  An object of class `formula` (or one that can be coerced to that
  class): a symbolic description of the model to be fitted. The details
  of model specification are given in 'Details'.

- ...:

  Additional `formula` objects to specify predictors of non-linear and
  distributional parameters. Formulas can either be named directly or
  contain names on their left-hand side. Alternatively, it is possible
  to fix parameters to certain values by passing numbers or character
  strings in which case arguments have to be named to provide the
  parameter names. See 'Details' for more information.

- flist:

  Optional list of formulas, which are treated in the same way as
  formulas passed via the `...` argument.

- family:

  Same argument as in
  [`brm`](https://paulbuerkner.com/brms/reference/brm.md). If `family`
  is specified in `brmsformula`, it will overwrite the value specified
  in other functions.

- autocor:

  An optional `formula` which contains autocorrelation terms as
  described in
  [`autocor-terms`](https://paulbuerkner.com/brms/reference/autocor-terms.md)
  or alternatively a
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md)
  object (deprecated). If `autocor` is specified in `brmsformula`, it
  will overwrite the value specified in other functions.

- nl:

  Logical; Indicates whether `formula` should be treated as specifying a
  non-linear model. By default, `formula` is treated as an ordinary
  linear model formula.

- loop:

  Logical; Only used in non-linear models. Indicates if the computation
  of the non-linear formula should be done inside (`TRUE`) or outside
  (`FALSE`) a loop over observations. Defaults to `TRUE`.

- center:

  Logical; Indicates if the population-level design matrix should be
  centered, which usually increases sampling efficiency. See the
  'Details' section for more information. Defaults to `TRUE` for
  distributional parameters and to `FALSE` for non-linear parameters.

- cmc:

  Logical; Indicates whether automatic cell-mean coding should be
  enabled when removing the intercept by adding `0` to the right-hand of
  model formulas. Defaults to `TRUE` to mirror the behavior of standard
  R formula parsing.

- sparse:

  Logical; indicates whether the population-level design matrices should
  be treated as sparse (defaults to `FALSE`). For design matrices with
  many zeros, this can considerably reduce required memory. Sampling
  speed is currently not improved or even slightly decreased.

- decomp:

  Optional name of the decomposition used for the population-level
  design matrix. Defaults to `NULL` that is no decomposition. Other
  options currently available are `"QR"` for the QR decomposition that
  helps in fitting models with highly correlated predictors.

- unused:

  An optional `formula` which contains variables that are unused in the
  model but should still be stored in the model's data frame. This can
  be useful, for example, if those variables are required for
  post-processing the model.

## Value

An object of class `brmsformula`, which is essentially a `list`
containing all model formulas as well as some additional information.

## Details

**General formula structure**

The `formula` argument accepts formulas of the following syntax:

`response | aterms ~ pterms + (gterms | group)`

The `pterms` part contains effects that are assumed to be the same
across observations. We call them 'population-level' or 'overall'
effects, or (adopting frequentist vocabulary) 'fixed' effects. The
optional `gterms` part may contain effects that are assumed to vary
across grouping variables specified in `group`. We call them
'group-level' or 'varying' effects, or (adopting frequentist vocabulary)
'random' effects, although the latter name is misleading in a Bayesian
context. For more details type `vignette("brms_overview")` and
`vignette("brms_multilevel")`.

**Group-level terms**

Multiple grouping factors each with multiple group-level effects are
possible. (Of course we can also run models without any group-level
effects.) Instead of `|` you may use `||` in grouping terms to prevent
correlations from being modeled. Equivalently, the `cor` argument of the
[`gr`](https://paulbuerkner.com/brms/reference/gr.md) function can be
used for this purpose, for example, `(1 + x || g)` is equivalent to
`(1 + x | gr(g, cor = FALSE))`.

It is also possible to model different group-level terms of the same
grouping factor as correlated (even across different formulas, e.g., in
non-linear models) by using `|<ID>|` instead of `|`. All group-level
terms sharing the same ID will be modeled as correlated. If, for
instance, one specifies the terms `(1+x|i|g)` and `(1+z|i|g)` somewhere
in the formulas passed to `brmsformula`, correlations between the
corresponding group-level effects will be estimated. In the above
example, `i` is not a variable in the data but just a symbol to indicate
correlations between multiple group-level terms. Equivalently, the `id`
argument of the [`gr`](https://paulbuerkner.com/brms/reference/gr.md)
function can be used as well, for example, `(1 + x | gr(g, id = "i"))`.

If levels of the grouping factor belong to different sub-populations, it
may be reasonable to assume a different covariance matrix for each of
the sub-populations. For instance, the variation within the treatment
group and within the control group in a randomized control trial might
differ. Suppose that `y` is the outcome, and `x` is the factor
indicating the treatment and control group. Then, we could estimate
different hyper-parameters of the varying effects (in this case a
varying intercept) for treatment and control group via
`y ~ x + (1 | gr(subject, by = x))`.

You can specify multi-membership terms using the
[`mm`](https://paulbuerkner.com/brms/reference/mm.md) function. For
instance, a multi-membership term with two members could be
`(1 | mm(g1, g2))`, where `g1` and `g2` specify the first and second
member, respectively. Moreover, if a covariate `x` varies across the
levels of the grouping-factors `g1` and `g2`, we can save the respective
covariate values in the variables `x1` and `x2` and then model the
varying effect as `(1 + mmc(x1, x2) | mm(g1, g2))`.

**Special predictor terms**

Flexible non-linear smooth terms can modeled using the
[`s`](https://paulbuerkner.com/brms/reference/s.md) and
[`t2`](https://paulbuerkner.com/brms/reference/s.md) functions in the
`pterms` part of the model formula. This allows to fit generalized
additive mixed models (GAMMs) with brms. The implementation is similar
to that used in the gamm4 package. For more details on this model class
see [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and
[`gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html).

Gaussian process terms can be fitted using the
[`gp`](https://paulbuerkner.com/brms/reference/gp.md) function in the
`pterms` part of the model formula. Similar to smooth terms, Gaussian
processes can be used to model complex non-linear relationships, for
instance temporal or spatial autocorrelation. However, they are
computationally demanding and are thus not recommended for very large
datasets or approximations need to be used.

The `pterms` and `gterms` parts may contain four non-standard effect
types namely monotonic, measurement error, missing value, and category
specific effects, which can be specified using terms of the form
`mo(predictor)`, `me(predictor, sd_predictor)`, `mi(predictor)`, and
`cs(<predictors>)`, respectively. Category specific effects can only be
estimated in ordinal models and are explained in more detail in the
package's main vignette (type `vignette("brms_overview")`). The other
three effect types are explained in the following.

A monotonic predictor must either be integer valued or an ordered
factor, which is the first difference to an ordinary continuous
predictor. More importantly, predictor categories (or integers) are not
assumed to be equidistant with respect to their effect on the response
variable. Instead, the distance between adjacent predictor categories
(or integers) is estimated from the data and may vary across categories.
This is realized by parameterizing as follows: One parameter takes care
of the direction and size of the effect similar to an ordinary
regression parameter, while an additional parameter vector estimates the
normalized distances between consecutive predictor categories. A main
application of monotonic effects are ordinal predictors that can this
way be modeled without (falsely) treating them as continuous or as
unordered categorical predictors. For more details and examples see
[`vignette("brms_monotonic")`](https://paulbuerkner.com/brms/articles/brms_monotonic.md).

Quite often, predictors are measured and as such naturally contain
measurement error. Although most researchers are well aware of this
problem, measurement error in predictors is ignored in most regression
analyses, possibly because only few packages allow for modeling it.
Notably, measurement error can be handled in structural equation models,
but many more general regression models (such as those featured by brms)
cannot be transferred to the SEM framework. In brms, effects of
noise-free predictors can be modeled using the `me` (for 'measurement
error') function. If, say, `y` is the response variable and `x` is a
measured predictor with known measurement error `sdx`, we can simply
include it on the right-hand side of the model formula via
`y ~ me(x, sdx)`. This can easily be extended to more general formulas.
If `x2` is another measured predictor with corresponding error `sdx2`
and `z` is a predictor without error (e.g., an experimental setting), we
can model all main effects and interactions of the three predictors in
the well known manner: `y ~ me(x, sdx) * me(x2, sdx2) * z`. The `me`
function is soft deprecated in favor of the more flexible and consistent
`mi` function (see below).

When a variable contains missing values, the corresponding rows will be
excluded from the data by default (row-wise exclusion). However, quite
often we want to keep these rows and instead estimate the missing
values. There are two approaches for this: (a) Impute missing values
before the model fitting for instance via multiple imputation (see
[`brm_multiple`](https://paulbuerkner.com/brms/reference/brm_multiple.md)
for a way to handle multiple imputed datasets). (b) Impute missing
values on the fly during model fitting. The latter approach is explained
in the following. Using a variable with missing values as predictors
requires two things, First, we need to specify that the predictor
contains missings that should to be imputed. If, say, `y` is the primary
response, `x` is a predictor with missings and `z` is a predictor
without missings, we go for `y ~ mi(x) + z`. Second, we need to model
`x` as an additional response with corresponding predictors and the
addition term [`mi()`](https://paulbuerkner.com/brms/reference/mi.md).
In our example, we could write `x | mi() ~ z`. Measurement error may be
included via the `sdy` argument, say, `x | mi(sdy = se) ~ z`. See
[`mi`](https://paulbuerkner.com/brms/reference/mi.md) for examples with
real data.

**Autocorrelation terms**

Autocorrelation terms can be directly specified inside the `pterms` part
as well. Details can be found in
[`autocor-terms`](https://paulbuerkner.com/brms/reference/autocor-terms.md).

**Additional response information**

Another special of the brms formula syntax is the optional `aterms`
part, which may contain multiple terms of the form `fun(<variable>)`
separated by `+` each providing special information on the response
variable. `fun` can be replaced with either `se`, `weights`, `subset`,
`cens`, `trunc`, `trials`, `cat`, `dec`, `rate`, `vreal`, or `vint`.
Their meanings are explained below (see also
[`addition-terms`](https://paulbuerkner.com/brms/reference/addition-terms.md)).

For families `gaussian`, `student` and `skew_normal`, it is possible to
specify standard errors of the observations, thus allowing to perform
meta-analysis. Suppose that the variable `yi` contains the effect sizes
from the studies and `sei` the corresponding standard errors. Then,
fixed and random effects meta-analyses can be conducted using the
formulas `yi | se(sei) ~ 1` and `yi | se(sei) ~ 1 + (1|study)`,
respectively, where `study` is a variable uniquely identifying every
study. If desired, meta-regression can be performed via
`yi | se(sei) ~ 1 + mod1 + mod2 + (1|study)` or  
`yi | se(sei) ~ 1 + mod1 + mod2 + (1 + mod1 + mod2|study)`, where `mod1`
and `mod2` represent moderator variables. By default, the standard
errors replace the parameter `sigma`. To model `sigma` in addition to
the known standard errors, set argument `sigma` in function `se` to
`TRUE`, for instance, `yi | se(sei, sigma = TRUE) ~ 1`.

For all families, weighted regression may be performed using `weights`
in the `aterms` part. Internally, this is implemented by multiplying the
log-posterior values of each observation by their corresponding weights.
Suppose that variable `wei` contains the weights and that `yi` is the
response variable. Then, formula `yi | weights(wei) ~ predictors`
implements a weighted regression.

For multivariate models, `subset` may be used in the `aterms` part, to
use different subsets of the data in different univariate models. For
instance, if `sub` is a logical variable and `y` is the response of one
of the univariate models, we may write `y | subset(sub) ~ predictors` so
that `y` is predicted only for those observations for which `sub`
evaluates to `TRUE`.

For log-linear models such as Poisson models, `rate` may be used in the
`aterms` part to specify the denominator of a response that is expressed
as a rate. The numerator is given by the actual response variable and
has a distribution according to the family as usual. In Poisson models,
using `rate(denom)` is equivalent to adding `offset(log(denom))` to the
linear predictor of the main parameter. For negative-binomial models,
this equivalence no longer holds and `rate(denom)` remains the
statistically correct approach.

With the exception of categorical and ordinal families, left, right, and
interval censoring can be modeled through
`y | cens(censored) ~ predictors`. The censoring variable (named
`censored` in this example) should contain the values `'left'`,
`'none'`, `'right'`, and `'interval'` (or equivalently `-1`, `0`, `1`,
and `2`) to indicate that the corresponding observation is left
censored, not censored, right censored, or interval censored. For
interval censored data, a second variable (let's call it `y2`) has to be
passed to `cens`. In this case, the formula has the structure
`y | cens(censored, y2) ~ predictors`. While the lower bounds are given
in `y`, the upper bounds are given in `y2` for interval censored data.
Intervals are assumed to be open on the left and closed on the right:
`(y, y2]`.

With the exception of categorical and ordinal families, the response
distribution can be truncated using the `trunc` function in the addition
part. If the response variable is truncated between, say, 0 and 100, we
can specify this via `yi | trunc(lb = 0, ub = 100) ~ predictors`.
Instead of numbers, variables in the data set can also be passed
allowing for varying truncation points across observations. Defining
only one of the two arguments in `trunc` leads to one-sided truncation.

For all continuous families, missing values in the responses can be
imputed within Stan by using the addition term `mi`. This is mostly
useful in combination with `mi` predictor terms as explained above under
'Special predictor terms'.

For families `binomial` and `zero_inflated_binomial`, addition should
contain a variable indicating the number of trials underlying each
observation. In `lme4` syntax, we may write for instance
`cbind(success, n - success)`, which is equivalent to
`success | trials(n)` in brms syntax. If the number of trials is
constant across all observations, say `10`, we may also write
`success | trials(10)`. **Please note that the
[`cbind()`](https://amices.org/mice/reference/cbind.html) syntax will
not work in brms in the expected way because this syntax is reserved for
other purposes.**

For all ordinal families, `aterms` may contain a term `thres(number)` to
specify the number thresholds (e.g, `thres(6)`), which should be equal
to the total number of response categories - 1. If not given, the number
of thresholds is calculated from the data. If different threshold
vectors should be used for different subsets of the data, the `gr`
argument can be used to provide the grouping variable (e.g,
`thres(6, gr = item)`, if `item` is the grouping variable). In this
case, the number of thresholds can also be a variable in the data with
different values per group.

A deprecated quasi alias of
[`thres()`](https://paulbuerkner.com/brms/reference/addition-terms.md)
is [`cat()`](https://paulbuerkner.com/brms/reference/addition-terms.md)
with which the total number of response categories (i.e., number of
thresholds + 1) can be specified.

In Wiener diffusion models (family `wiener`) the addition term `dec` is
mandatory to specify the (vector of) binary decisions corresponding to
the reaction times. Non-zero values will be treated as a response on the
upper boundary of the diffusion process and zeros will be treated as a
response on the lower boundary. Alternatively, the variable passed to
`dec` might also be a character vector consisting of `'lower'` and
`'upper'`.

All families support the `index` addition term to uniquely identify each
observation of the corresponding response variable. Currently, `index`
is primarily useful in combination with the `subset` addition and
[`mi`](https://paulbuerkner.com/brms/reference/mi.md) terms.

For custom families, it is possible to pass an arbitrary number of real
and integer vectors via the addition terms `vreal` and `vint`,
respectively. An example is provided in
[`vignette('brms_customfamilies')`](https://paulbuerkner.com/brms/articles/brms_customfamilies.md).
To pass multiple vectors of the same data type, provide them separated
by commas inside a single `vreal` or `vint` statement.

Multiple addition terms of different types may be specified at the same
time using the `+` operator. For example, the formula
`formula = yi | se(sei) + cens(censored) ~ 1` implies a censored
meta-analytic model.

The addition argument `disp` (short for dispersion) has been removed in
version 2.0. You may instead use the distributional regression approach
by specifying `sigma ~ 1 + offset(log(xdisp))` or
`shape ~ 1 + offset(log(xdisp))`, where `xdisp` is the variable being
previously passed to `disp`.

**Parameterization of the population-level intercept**

By default, the population-level intercept (if incorporated) is
estimated separately and not as part of population-level parameter
vector `b` As a result, priors on the intercept also have to be
specified separately. Furthermore, to increase sampling efficiency, the
population-level design matrix `X` is centered around its column means
`X_means` if the intercept is incorporated. This leads to a temporary
bias in the intercept equal to `<X_means, b>`, where `<,>` is the scalar
product. The bias is corrected after fitting the model, but be aware
that you are effectively defining a prior on the intercept of the
centered design matrix not on the real intercept. You can turn off this
special handling of the intercept by setting argument `center` to
`FALSE`. For more details on setting priors on population-level
intercepts, see
[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md).

This behavior can be avoided by using the reserved (and internally
generated) variable `Intercept`. Instead of `y ~ x`, you may write
`y ~ 0 + Intercept + x`. This way, priors can be defined on the real
intercept, directly. In addition, the intercept is just treated as an
ordinary population-level effect and thus priors defined on `b` will
also apply to it. Note that this parameterization may be less efficient
than the default parameterization discussed above.

**Formula syntax for non-linear models**

In brms, it is possible to specify non-linear models of arbitrary
complexity. The non-linear model can just be specified within the
`formula` argument. Suppose, that we want to predict the response `y`
through the predictor `x`, where `x` is linked to `y` through
`y = alpha - beta * lambda^x`, with parameters `alpha`, `beta`, and
`lambda`. This is certainly a non-linear model being defined via
`formula = y ~ alpha - beta * lambda^x` (addition arguments can be added
in the same way as for ordinary formulas). To tell brms that this is a
non-linear model, we set argument `nl` to `TRUE`. Now we have to specify
a model for each of the non-linear parameters. Let's say we just want to
estimate those three parameters with no further covariates or random
effects. Then we can pass `alpha + beta + lambda ~ 1` or equivalently
(and more flexible) `alpha ~ 1, beta ~ 1, lambda ~ 1` to the `...`
argument. This can, of course, be extended. If we have another predictor
`z` and observations nested within the grouping factor `g`, we may write
for instance `alpha ~ 1, beta ~ 1 + z + (1|g), lambda ~ 1`. The formula
syntax described above applies here as well. In this example, we are
using `z` and `g` only for the prediction of `beta`, but we might also
use them for the other non-linear parameters (provided that the
resulting model is still scientifically reasonable).

By default, non-linear covariates are treated as real vectors in Stan.
However, if the data of the covariates is of type \`integer\` in R
(which can be enforced by the \`as.integer\` function), the Stan type
will be changed to an integer array. That way, covariates can also be
used for indexing purposes in Stan.

Parameters of non-linear models may not be uniquely identified or the
posterior inference may have convergence issues. For this reason it is
mandatory to specify priors on the non-linear parameters. For
instructions on how to do that, see
[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md). For
some examples of non-linear models, see
[`vignette("brms_nonlinear")`](https://paulbuerkner.com/brms/articles/brms_nonlinear.md).

**Formula syntax for predicting distributional parameters**

It is also possible to predict parameters of the response distribution
such as the residual standard deviation `sigma` in gaussian models or
the hurdle probability `hu` in hurdle models. The syntax closely
resembles that of a non-linear parameter, for instance
`sigma ~ x + s(z) + (1+x|g)`. For some examples of distributional
models, see
[`vignette("brms_distreg")`](https://paulbuerkner.com/brms/articles/brms_distreg.md).

Parameter `mu` exists for every family and can be used as an alternative
to specifying terms in `formula`. If both `mu` and `formula` are given,
the right-hand side of `formula` is ignored. Accordingly, specifying
terms on the right-hand side of both `formula` and `mu` at the same time
is deprecated. In future versions, `formula` might be updated by `mu`.

The following are distributional parameters of specific families (all
other parameters are treated as non-linear parameters): `sigma`
(residual standard deviation or scale of the `gaussian`, `student`,
`skew_normal`, `lognormal` `exgaussian`, and `asym_laplace` families);
`shape` (shape parameter of the `Gamma`, `weibull`, `negbinomial`, and
related zero-inflated / hurdle families); `nu` (degrees of freedom
parameter of the `student` and `frechet` families); `phi` (precision
parameter of the `beta`, `zero_inflated_beta`, and `xbeta` families);
`kappa` (precision parameter of the `von_mises` family); `beta` (mean
parameter of the exponential component of the `exgaussian` family);
`quantile` (quantile parameter of the `asym_laplace` family); `zi`
(zero-inflation probability); `hu` (hurdle probability); `zoi`
(zero-one-inflation probability); `coi` (conditional one-inflation
probability); `disc` (discrimination) for ordinal models; `bs`, `ndt`,
and `bias` (boundary separation, non-decision time, and initial bias of
the `wiener` diffusion model). By default, distributional parameters are
modeled on the log scale if they can be positive only or on the logit
scale if the can only be within the unit interval.

Alternatively, one may fix distributional parameters to certain values.
However, this is mainly useful when a model becomes too complicated and
posterior inference has convergence issues. We thus suggest to be
generally careful when making use of this option. The `quantile`
parameter of the `asym_laplace` distribution is a good example where it
is useful. By fixing `quantile`, one can perform quantile regression for
the specified quantile. For instance, `quantile = 0.25` allows
predicting the 25%-quantile. Furthermore, the `bias` parameter in
drift-diffusion models, is assumed to be `0.5` (i.e. no bias) in many
applications. To achieve this, simply write `bias = 0.5`. Other possible
applications are the Cauchy distribution as a special case of the
Student-t distribution with `nu = 1`, or the geometric distribution as a
special case of the negative binomial distribution with `shape = 1`.
Furthermore, the parameter `disc` ('discrimination') in ordinal models
is fixed to `1` by default and not estimated, but may be modeled as any
other distributional parameter if desired (see examples). For reasons of
identification, `'disc'` can only be positive, which is achieved by
applying the log-link.

In categorical models, distributional parameters do not have fixed
names. Instead, they are named after the response categories (excluding
the first one, which serves as the reference category), with the prefix
`'mu'`. If, for instance, categories are named `cat1`, `cat2`, and
`cat3`, the distributional parameters will be named `mucat2` and
`mucat3`.

Some distributional parameters currently supported by `brmsformula` have
to be positive (a negative standard deviation or precision parameter
does not make any sense) or are bounded between 0 and 1 (for
zero-inflated / hurdle probabilities, quantiles, or the initial bias
parameter of drift-diffusion models). However, linear predictors can be
positive or negative, and thus the log link (for positive parameters) or
logit link (for probability parameters) are used by default to ensure
that distributional parameters are within their valid intervals. This
implies that, by default, effects for such distributional parameters are
estimated on the log / logit scale and one has to apply the inverse link
function to get to the effects on the original scale. Alternatively, it
is possible to use the identity link to predict parameters on their
original scale, directly. However, this is much more likely to lead to
problems in the model fitting, if the parameter actually has a
restricted range.

See also
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md)
for an overview of valid link functions.

**Formula syntax for mixture models**

The specification of mixture models closely resembles that of
non-mixture models. If not specified otherwise (see below), all mean
parameters of the mixture components are predicted using the right-hand
side of `formula`. All types of predictor terms allowed in non-mixture
models are allowed in mixture models as well.

Distributional parameters of mixture distributions have the same name as
those of the corresponding ordinary distributions, but with a number at
the end to indicate the mixture component. For instance, if you use
family `mixture(gaussian, gaussian)`, the distributional parameters are
`sigma1` and `sigma2`. Distributional parameters of the same class can
be fixed to the same value. For the above example, we could write
`sigma2 = "sigma1"` to make sure that both components have the same
residual standard deviation, which is in turn estimated from the data.

In addition, there are two types of special distributional parameters.
The first are named `mu<ID>`, that allow for modeling different
predictors for the mean parameters of different mixture components. For
instance, if you want to predict the mean of the first component using
predictor `x` and the mean of the second component using predictor `z`,
you can write `mu1 ~ x` as well as `mu2 ~ z`. The second are named
`theta<ID>`, which constitute the mixing proportions. If the mixing
proportions are fixed to certain values, they are internally normalized
to form a probability vector. If one seeks to predict the mixing
proportions, all but one of the them has to be predicted, while the
remaining one is used as the reference category to identify the model.
The so-called 'softmax' transformation is applied on the linear
predictor terms to form a probability vector.

For more information on mixture models, see the documentation of
[`mixture`](https://paulbuerkner.com/brms/reference/mixture.md).

**Formula syntax for multivariate models**

Multivariate models may be specified using `mvbind` notation or with
help of the
[`mvbf`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)
function. Suppose that `y1` and `y2` are response variables and `x` is a
predictor. Then `mvbind(y1, y2) ~ x` specifies a multivariate model. The
effects of all terms specified at the RHS of the formula are assumed to
vary across response variables. For instance, two parameters will be
estimated for `x`, one for the effect on `y1` and another for the effect
on `y2`. This is also true for group-level effects. When writing, for
instance, `mvbind(y1, y2) ~ x + (1+x|g)`, group-level effects will be
estimated separately for each response. To model these effects as
correlated across responses, use the ID syntax (see above). For the
present example, this would look as follows:
`mvbind(y1, y2) ~ x + (1+x|2|g)`. Of course, you could also use any
value other than `2` as ID.

It is also possible to specify different formulas for different
responses. If, for instance, `y1` should be predicted by `x` and `y2`
should be predicted by `z`, we could write `mvbf(y1 ~ x, y2 ~ z)`.
Alternatively, multiple `brmsformula` objects can be added to specify a
joint multivariate model (see 'Examples').

## See also

[`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md),
[`brmsformula-helpers`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.md)

## Examples

``` r
# multilevel model with smoothing terms
brmsformula(y ~ x1*x2 + s(z) + (1+x1|1) + (1|g2))
#> y ~ x1 * x2 + s(z) + (1 + x1 | 1) + (1 | g2) 

# additionally predict 'sigma'
brmsformula(y ~ x1*x2 + s(z) + (1+x1|1) + (1|g2),
            sigma ~ x1 + (1|g2))
#> y ~ x1 * x2 + s(z) + (1 + x1 | 1) + (1 | g2) 
#> sigma ~ x1 + (1 | g2)

# use the shorter alias 'bf'
(formula1 <- brmsformula(y ~ x + (x|g)))
#> y ~ x + (x | g) 
(formula2 <- bf(y ~ x + (x|g)))
#> y ~ x + (x | g) 
# will be TRUE
identical(formula1, formula2)
#> [1] TRUE

# incorporate censoring
bf(y | cens(censor_variable) ~ predictors)
#> y | cens(censor_variable) ~ predictors 

# define a simple non-linear model
bf(y ~ a1 - a2^x, a1 + a2 ~ 1, nl = TRUE)
#> y ~ a1 - a2^x 
#> a1 ~ 1
#> a2 ~ 1

# predict a1 and a2 differently
bf(y ~ a1 - a2^x, a1 ~ 1, a2 ~ x + (x|g), nl = TRUE)
#> y ~ a1 - a2^x 
#> a1 ~ 1
#> a2 ~ x + (x | g)

# correlated group-level effects across parameters
bf(y ~ a1 - a2^x, a1 ~ 1 + (1 |2| g), a2 ~ x + (x |2| g), nl = TRUE)
#> y ~ a1 - a2^x 
#> a1 ~ 1 + (1 | 2 | g)
#> a2 ~ x + (x | 2 | g)
# alternative but equivalent way to specify the above model
bf(y ~ a1 - a2^x, a1 ~ 1 + (1 | gr(g, id = 2)),
   a2 ~ x + (x | gr(g, id = 2)), nl = TRUE)
#> y ~ a1 - a2^x 
#> a1 ~ 1 + (1 | gr(g, id = 2))
#> a2 ~ x + (x | gr(g, id = 2))

# define a multivariate model
bf(mvbind(y1, y2) ~ x * z + (1|g))
#> y1 ~ x * z + (1 | g) 
#> y2 ~ x * z + (1 | g) 

# define a zero-inflated model
# also predicting the zero-inflation part
bf(y ~ x * z + (1+x|ID1|g), zi ~ x + (1|ID1|g))
#> y ~ x * z + (1 + x | ID1 | g) 
#> zi ~ x + (1 | ID1 | g)

# specify a predictor as monotonic
bf(y ~ mo(x) + more_predictors)
#> y ~ mo(x) + more_predictors 

# for ordinal models only
# specify a predictor as category specific
bf(y ~ cs(x) + more_predictors)
#> y ~ cs(x) + more_predictors 
# add a category specific group-level intercept
bf(y ~ cs(x) + (cs(1)|g))
#> y ~ cs(x) + (cs(1) | g) 
# specify parameter 'disc'
bf(y ~ person + item, disc ~ item)
#> y ~ person + item 
#> disc ~ item

# specify variables containing measurement error
bf(y ~ me(x, sdx))
#> y ~ me(x, sdx) 

# specify predictors on all parameters of the wiener diffusion model
# the main formula models the drift rate 'delta'
bf(rt | dec(decision) ~ x, bs ~ x, ndt ~ x, bias ~ x)
#> rt | dec(decision) ~ x 
#> bs ~ x
#> ndt ~ x
#> bias ~ x

# fix the bias parameter to 0.5
bf(rt | dec(decision) ~ x, bias = 0.5)
#> rt | dec(decision) ~ x 
#> bias = 0.5

# specify different predictors for different mixture components
mix <- mixture(gaussian, gaussian)
#> Setting order = 'mu' for mixtures of the same family.
bf(y ~ 1, mu1 ~ x, mu2 ~ z, family = mix)
#> y ~ 1 
#> mu1 ~ x
#> mu2 ~ z

# fix both residual standard deviations to the same value
bf(y ~ x, sigma2 = "sigma1", family = mix)
#> y ~ x 
#> sigma2 = sigma1

# use the '+' operator to specify models
bf(y ~ 1) +
  nlf(sigma ~ a * exp(b * x), a ~ x) +
  lf(b ~ z + (1|g), dpar = "sigma") +
  gaussian()
#> Warning: Arguments '...' and 'flist' in nlf() will be reworked at some point. Please avoid using them if possible.
#> Warning: Argument 'dpar' is no longer necessary and ignored.
#> y ~ 1 
#> sigma ~ a * exp(b * x)
#> a ~ x
#> b ~ z + (1 | g)

# specify a multivariate model using the '+' operator
bf(y1 ~ x + (1|g)) +
  gaussian() + cor_ar(~1|g) +
  bf(y2 ~ z) + poisson()
#> Warning: Using 'cor_brms' objects for 'autocor' is deprecated. Please see ?cor_brms for details.
#> y1 ~ x + (1 | g) 
#> autocor ~ arma(time = NA, gr = g, p = 1, q = 0, cov = FALSE)
#> y2 ~ z 

# specify correlated residuals of a gaussian and a poisson model
form1 <- bf(y1 ~ 1 + x + (1|c|obs), sigma = 1) + gaussian()
form2 <- bf(y2 ~ 1 + x + (1|c|obs)) + poisson()

# model missing values in predictors
bf(bmi ~ age * mi(chl)) +
  bf(chl | mi() ~ age) +
  set_rescor(FALSE)
#> bmi ~ age * mi(chl) 
#> chl | mi() ~ age 

# model sigma as a function of the mean
bf(y ~ eta, nl = TRUE) +
  lf(eta ~ 1 + x) +
  nlf(sigma ~ tau * sqrt(eta)) +
  lf(tau ~ 1)
#> y ~ eta 
#> eta ~ 1 + x
#> sigma ~ tau * sqrt(eta)
#> tau ~ 1
```
