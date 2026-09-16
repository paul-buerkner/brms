# Physical sum-to-zero implementation

The public option `gr(..., s2z = TRUE)` changes sampling coordinates while
preserving the conventional population and group-effect model. The supported
formulas and priors are documented in `?gr`. This note describes the internal
contract and the invariants that the implementation must preserve.

## Coordinates and reconstruction

For each group-effect block, write the conventional effects as
`r = r_s2z + 1 * m`, where each column of `r_s2z` sums to zero over observed
grouping levels and `m` contains the omitted coefficient means. Each varying
design column must match a population design column, so its omitted mean can
be absorbed into a population coordinate. Writing the combined map as `A`,
the likelihood uses `theta = b + A * m` and `r_s2z`.

The population design's ordinary intercept centering is part of this map.
Consequently, identifying an affected population coefficient requires the
design map, not just a check for a matching coefficient label. Multiple
blocks in one predictor share a population coefficient system and their
omitted means must be handled jointly.

Conventional `b` and `r` are reconstructed in generated quantities. Public
`fixef()`, `ranef()`, `coef()`, and prediction operate on those conventional
draws. Reported group effects therefore need not sum to zero. Public
prediction, SAR, GP, and autocorrelation implementations require no S2Z
algorithm. Custom Stan code cannot use the reconstructed parameters before
generated quantities.

## Plan and code generation

`R/re-s2z.R` owns validation and the predictor-local plan. `re_s2z_plan()`
returns `NULL` for an ordinary predictor. An S2Z plan records population
coordinate names and original indices, centered-design metadata, active and
inactive coordinate maps, and the participating blocks. A block records its
group frame, matching population indices, active-coordinate indices, and
diagnostic context.

Structural plans can be cached in `frame$re_s2z_plan`. Prior-bearing plans
must be resolved from the effective prior table supplied by the current
call: prior information is never stored in the structural cache. The resolved
`prior` list is restricted to active coordinates and named by their
population coefficient names. Fixed-only coordinates retain original
population indices for ordinary prior expressions, including vector-valued
arguments.

`R/stan-re-s2z.R` owns S2Z code generation. `stan_re_s2z_systems()` obtains
the plan and chooses the scalar, independent, dense, explicit-mean, Matheron,
or joint path. `stan_re_s2z_block()` and `stan_re_s2z_system()` provide the
block and system dispatch boundaries. Singleton emission retains its
established ordering, including generated-quantities RNG calls.

`stan_fe_s2z()` handles population coordinates. `stan_re_s2z_coef()` supplies
the likelihood coefficient expression to both ordinary matrix multiplication
and GLM primitives. `stan_re_s2z_public_fe_def()` and
`stan_re_s2z_public_fe_comp()` own the declarations and reconstruction of
public population effects.

`stan_re_s2z_public_re_comp()` reconstructs public group effects from the
physical deviations and the kernel's scalar or vector mean. It preserves the
matrix/column ordering required by correlated effects. Internal-output names
are kept alongside these generators; ordinary `save_pars()` handling delegates
to `re_s2z_internal_fe()` and `re_s2z_internal_re()`.

`re_s2z_plan_infos()` is a compatibility adapter for the existing specialized
generators. It produces their full population-coordinate prior list, with
flat placeholders for fixed-only coordinates, and their repeated metadata.
This duplication is transient; the plan itself remains the source of truth.
New generator interfaces should consume plan fields directly rather than
extend the compatibility representation.

## Priors and specialized paths

Only population coordinates affected by an omitted-mean map enter the S2Z
prior calculation. Their supported flat, normal, Student-t, Cauchy, and
logistic priors retain their exact densities. Group-level scale and
correlation priors retain their conventional interpretation; Student-t `sd`
parameters are scales, not marginal standard deviations.

Conditional Gaussian systems integrate omitted means analytically and draw
them from their conditional distribution during reconstruction. Independent
and scalar cases use the corresponding simpler algebra. The Gaussian
Matheron path uses the same target model with a specialized computation.
Student-t effects and population priors retain the scale-mixture variables
required by their exact densities.

Logistic-active systems explicitly sample standardized omitted means and
include the exact logistic density and change-of-variables contribution.
This supports default logistic intercept priors in auxiliary probability
predictors without a Gaussian approximation. Fixed-only coordinates retain
ordinary declarations, priors, and vector argument indexing.

## Validation boundaries

Validate structure and matching designs when constructing frames, and
validate prior restrictions after ordinary and special priors have been
resolved. Public `validate_prior()`, `standata()`, and `stancode()` must agree
about supported models. Unsupported structures must fail before compilation
with enough predictor, group, and coefficient context to correct the model.

Tests should establish the zero-sum geometry, transformed-density and
reconstruction equivalence, preservation of original prior indices, and
unchanged ordinary-model code generation. Native model comparisons cover
posterior quantities and predictive behavior. Specialized performance paths
must preserve target values, gradients, and the public reconstruction
contract. Generic prior-generation and prediction-metadata bugs are separate
prerequisites with ordinary-model regression tests.

`tests/local/tests.s2z-compatibility.R` compares separately installed revisions
in fresh R processes. Use `conventional` mode against the upstream revision
plus the explicit prerequisites, and `s2z` mode against the previous foundation.
It requires exact Stan source and data equality; the S2Z comparison permits
only deletion of the checked, unused historical `chol2inv_brms` definition.
It also compares diagnostics for the unsupported sparse and QR S2Z designs.
