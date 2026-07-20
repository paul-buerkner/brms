# PR summary: extended `posterior_predict()` outputs

Branch: `posterior_pit`  
Related PR: [paul-buerkner/brms#1857](https://github.com/paul-buerkner/brms/pull/1857)  

## Detailed overviews (upstream PR comments)

Family-by-family support tables live in two comments on #1857:

- [`posterior_predict()` `output` support](https://github.com/paul-buerkner/brms/pull/1857#issuecomment-4325781243) — which methods support `random` / `probability` / `density` / `pit` / `quantile`
- [Distribution helpers in `R/distributions.R`](https://github.com/paul-buerkner/brms/pull/1857#issuecomment-4325820624) — which families have `d` / `p` / `q` / `r` (including newly added helpers)

## What

Extends `posterior_predict()` beyond random draws so it can return CDF/PIT values, densities/PMFs, and quantiles for many response families, with matching distribution helpers (`p`/`d`/`q`/`r`) where needed.

## Changes

### API (`posterior_predict()` / `predict()`)

New arguments (default preserves old behavior: `output = "random"`):

| Argument | Role |
|----------|------|
| `output` | `"random"`, `"probability"`, `"pit"`, `"density"`, or `"quantile"` |
| `q` | query point for probability / PIT / density (defaults to observed `Y`) |
| `p` | probability for quantiles |
| `lower.tail`, `log.p`, `log` | CDF / quantile / density options |

Shared backends: `predict_continuous_helper()`, `predict_discrete_helper()`, plus validation that rejects unsupported `output` × family combinations.

### Distributions (`R/distributions.R`)

New or extended helpers used by the predictive methods, including quantile functions for inverse Gaussian, ex-Gaussian, von Mises, COM-Poisson, beta-binomial, zero-inflated / hurdle families, and `d`/`p`/`q` for zero-one-inflated beta; plus ordinal / categorical / hurdle-cumulative support.

### Tests and docs

- Spec-driven tests (`helper-distributions.R`, `helper-posterior-predict.R`) covering supported modes, truncation where applicable, and relations between `d`/`p`/`q`/`r`
- Vignettes: `brms_posterior_predict.Rmd`, `brms_distribution_output_cheatsheet.Rmd`

## Caveats

1. **Not all families support all modes.** Multivariate / time / spatial / residual-correlation Gaussians and Students, Wiener, multinomial / Dirichlet(-multinomial) / logistic-normal, and `mixture` accept only `output = "random"`. `cox` still has no predictive draws. Custom families default to `random` unless they implement other modes.

2. **Discrete PIT is randomized.** Continuous `probability` and `pit` coincide; discrete `pit` uses randomized PIT to avoid point-mass spikes.

3. **Discrete truncation still uses rejection sampling** (`ntrys`); pathological bounds can be slow or approximate (existing brms behavior). Hurdle families now respect truncation for `output = "random"` (continuous via inverse CDF; discrete via rejection). Some zero-inflated families still overlay exact zeros after truncated base sampling, so draws of 0 can fall outside `lb > 0` bounds (pre-existing pattern).

4. **Some quantiles are numerical.** Families without closed-form inverses (e.g. inverse Gaussian, ex-Gaussian, von Mises, COM-Poisson, several zero-inflated / hurdle CDFs) invert the CDF numerically (`tol` where exposed).

5. **Optional Suggests.** Full coverage of `skew_normal` / `xbeta` in the cheatsheet needs `mnormt` / `betareg`.


