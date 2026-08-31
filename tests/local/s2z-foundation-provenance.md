# PLAN-01 S2Z foundation provenance

This manifest records how the combined development branch preserved at
`archive/s2z-combined-20260825` was reduced to the physical/shared-scale
foundation on `feature/s2z-foundation`. The original extraction base was
`29ccd2384776b2bfc13bf749659528911743a92b`; the resulting series was then
rebased onto `origin/master` at
`4f91054353436f5814cd4ccf6da61f14b2d0e45d`.

| Combined commit | Foundation treatment |
|---|---|
| `112e2fce` Add marginalized sum-to-zero group effects | Replayed as `562a969f`. |
| `fd3b08b3` Specialize scalar sum-to-zero effects | Replayed as `a7449e20`. |
| `1835a8a5` Document varying-effect prior support | Replayed as `003ea1e4`; this describes coefficient-specific scales shared across levels. |
| `4b36a563` Specialize independent sum-to-zero effects | Replayed as `c1afad55`. |
| `601e792c` Clarify physical sum-to-zero terminology | Replayed as `2a99dc40`. |
| `5f692168` Add group-varying S2Z scales | Deferred to PLAN-04; no public API or state retained. |
| `2a4eeb11` Add partial centering for sum-to-zero effects | Deferred to PLAN-03; no `center`, `rho_s2z`, or automatic-centering state retained. |
| `ca59f007` Support multiple sum-to-zero grouping blocks | Replayed as `48169204` after removing dependencies on deferred modes. |
| `4c95250a` Add fast Matheron path for Gaussian S2Z blocks | Replayed as `e1381595` after removing dependencies on deferred modes. |
| `50ec7242` Allow priors on realized group-level scales | Deferred to PLAN-04; no realized-scale API or level-specific scale state retained. |
| `ddf913c5` Reduce Stan allocations for S2Z effects | PR-1-safe allocation changes were rewritten onto the reduced physical paths. |
| `ce57b03b` Vectorize S2Z scale priors | Shared-scale vectorization was rewritten where applicable; varying-level-scale hunks were deferred to PLAN-04. |
| `9d76772e` Vectorize additional S2Z operations | PR-1-safe vectorization was rewritten; its varying-scale fixture was not retained. |
| `7aed341b` Add automatic Fisher centering for S2Z effects | Deferred to PLAN-03; no Fisher-centering API or generated state retained. |

PLAN-01 also adds work that postdates the combined branch: exact logistic
population priors, phased structure/design/prior validation, contextual
diagnostics, active versus fixed-only population coordinates, a semantic
capability matrix, and posterior-equivalence fixtures.

The mechanical forbidden-feature audit lives in
`tests/testthat/tests.s2z-capabilities.R`. The combined branch and its pushed
copy (`feature/s2z-group-effects`) remain unchanged and recoverable.
