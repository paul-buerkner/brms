# local test for equivalence between asym_laplace and skew_double_exponential (currently gain ~ 25% speed)
test_that("new skew_double_exponential is equivalent to asym_laplace", suppressWarnings({
  data(Mammals, package='quantreg')
  fit_asym <- brm(bf(weight~speed, quantile=0.8), family=asym_laplace(), iter=100000, warmup=500, refresh=0, cores=4, chains=4, seed=1, data=Mammals, backend='cmdstanr')
  fit_skew <- brm(bf(weight~speed, tau=0.8), family=skew_double_exponential(),  iter=100000, warmup=500, refresh=0, cores=4, chains=4, seed=1, data=Mammals, backend='cmdstanr')
  
  delta <- fixef(fit_asym)[,1] - fixef(fit_skew)[,1]
  expect_all_true(delta/fixef(fit_asym)[,1] < 1e-2)
}))
