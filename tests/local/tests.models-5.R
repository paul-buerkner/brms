source("setup_tests_local.R")

old_bayes_R2 <- function(fit) {
  y <- brms:::get_y(fit, warn = TRUE)
  ypred <- posterior_epred(fit, sort = TRUE)
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- matrixStats::rowVars(ypred)
  var_e <- matrixStats::rowVars(e)
  var_ypred / (var_ypred + var_e)
}

test_that("bayes_R2 prior-only draws match model-based formula", {
  set.seed(2801)
  dat <- data.frame(
    x = rnorm(40),
    y = rnorm(40)
  )
  fit <- brm(
    y ~ x, data = dat, family = gaussian(),
    prior = c(
      prior(normal(0, 0.5), class = b),
      prior(normal(0, 0.5), class = Intercept),
      prior(exponential(1), class = sigma)
    ),
    sample_prior = "only",
    chains = 1, iter = 600, warmup = 200,
    seed = 2802, refresh = 0
  )
  ypred <- posterior_epred(fit, sort = TRUE)
  sigma <- posterior_epred(fit, dpar = "sigma", sort = TRUE)
  var_mupred <- matrixStats::rowVars(ypred)
  var_res <- matrixStats::rowMeans2(sigma^2)
  expected <- var_mupred / (var_mupred + var_res)
  expect_equal(bayes_R2(fit, summary = FALSE)[, 1], expected)
})

test_that("bayes_R2 posterior difference shrinks with larger n", {
  make_fit <- function(n, seed) {
    set.seed(seed)
    dat <- data.frame(x = rnorm(n))
    dat$y <- 1 + 1.5 * dat$x + rnorm(n, sd = 1)
    brm(
      y ~ x, data = dat, family = gaussian(),
      chains = 1, iter = 700, warmup = 300,
      seed = seed, refresh = 0
    )
  }
  fit_small <- make_fit(8, 2803)
  fit_large <- make_fit(250, 2804)

  model_small <- bayes_R2(fit_small, summary = FALSE)[, 1]
  model_large <- bayes_R2(fit_large, summary = FALSE)[, 1]
  residual_small <- old_bayes_R2(fit_small)
  residual_large <- old_bayes_R2(fit_large)

  diff_small <- abs(stats::median(model_small) - stats::median(residual_small))
  diff_large <- abs(stats::median(model_large) - stats::median(residual_large))
  expect_gt(diff_small, diff_large)
})

test_that("global shrinkage priors work correctly", suppressWarnings({
  set.seed(2563)
  dat <- epilepsy
  dat$x1 <- sample(1:4, nrow(dat), TRUE)
  dat$x2 <- sample(1:10, nrow(dat), TRUE)
  dat$x3 <- sample(1:10, nrow(dat), TRUE)

  bprior <- prior(R2D2(main = TRUE), class = b) +
    prior(R2D2(), class = sd) +
    prior(R2D2(), class = sds) +
    prior(R2D2(), class = sdgp) +
    prior(R2D2(), class = ar) +
    prior(R2D2(), class = ma)
  bform <- bf(count ~ Trt * Base + Age + mo(x1) + (1|patient) +
                gp(x2) + s(x3) + arma(p = 2, q = 2, gr = patient))
  fit <- brm(bform, data = dat, prior = bprior, cores = 4,
             control = list(adapt_delta = 0.95), seed = 8892)

  classes <- c("sdb", "sdbsp", "sdbs", "sdar", "sdma")
  for (cl in classes) {
    expect_true(any(grepl(paste0("^", cl), variables(fit))))
  }
  expect_range(SW(waic(fit))$estimates[3, 1], 1580, 1660)
}))

test_that("Non-linear matrix covariates work correctly", suppressWarnings({
  set.seed(2134)
  N <- 100
  dat <- data.frame(y=rnorm(N))
  dat$X <- matrix(rnorm(N*2), N, 2)

  nlfun_stan <- "
    real nlfun(real a, real b, real c, row_vector X) {
       return a + b * X[1] + c * X[2];
    }
  "
  nlstanvar <- stanvar(scode = nlfun_stan, block = "functions")

  # fit the model
  bform <- bf(y~nlfun(a, b, c, X), a~1, b~1, c~1, nl = TRUE)
  fit <- brm(bform, dat, stanvars = nlstanvar)
  b <- as.matrix(fit, variable = "b_b_Intercept")
  # fit benchmark model that should yield the same results up to MCMC error
  fit2 <- brm(y~X, dat)
  b2 <- as.matrix(fit2, variable = "b_X1")
  expect_range(mean(b) - mean(b2), -0.1, 0.1)

  # test that post processing works too
  nlfun <- function(a, b, c, X) {
    a + b * X[, , 1] + c * X[, , 2]
  }
  expect_equal(dim(posterior_epred(fit)), c(ndraws(fit), nobs(fit)))
  expect_equal(dim(posterior_predict(fit)), c(ndraws(fit), nobs(fit)))

  lcomp <- loo(fit, fit2)$diffs
  expect_range(lcomp[2, 1], -1, 1)
}))

test_that("Addition argument 'subset' works correctly", suppressWarnings({
  set.seed(12454)
  data("BTdata", package = "MCMCglmm")
  BTdata$sub1 <- sample(0:1, nrow(BTdata), replace = TRUE)
  BTdata$sub2 <- sample(0:1, nrow(BTdata), replace = TRUE)

  bform <- bf(tarsus | subset(sub1) ~ sex + (1|p|fosternest) + (1|q|dam)) +
    bf(back | subset(sub2) ~ sex + (1|p|fosternest) + (1|q|dam)) +
    set_rescor(FALSE)
  fit <- brm(bform, BTdata, refresh = 0)
  print(summary(fit))
  expect_error(predict(fit), "'resp' must be a single variable name")
  pred <- predict(fit, resp = "tarsus")
  expect_equal(nrow(pred), sum(BTdata$sub1))
  pred <- fitted(fit, resp = "back")
  expect_equal(nrow(pred), sum(BTdata$sub2))
  waic <- waic(fit, resp = "back")
  expect_range(waic$estimates[3, 1], 1100, 1200)
  ce <- conditional_effects(fit)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])
  expect_equal(nobs(fit, resp = "tarsus"), sum(BTdata$sub1))
}))

test_that("Cox models work correctly", suppressWarnings({
  set.seed(12345)
  covs <- data.frame(id  = 1:200, trt = stats::rbinom(200, 1L, 0.5))
  d1 <- simsurv::simsurv(lambdas = 0.1, gammas  = 1.5, betas = c(trt = -0.5),
                         x = covs, maxt  = 5)
  d1 <- merge(d1, covs)
  d1$g <- sample(c("a", "b"), nrow(d1), TRUE)

  fit1 <- brm(eventtime | cens(1 - status) ~ 1 + trt,
              data = d1, family = brmsfamily("cox"), refresh = 0)
  print(summary(fit1))
  expect_range(posterior_summary(fit1)["b_trt", "Estimate"], -0.70, -0.30)
  expect_range(waic(fit1)$estimates[3, 1], 620, 670)

  fit2 <- brm(eventtime | cens(1 - status) + bhaz(gr = g) ~ 1 + trt,
              data = d1, family = brmsfamily("cox"), refresh = 0)
  print(summary(fit2))
  expect_true("sbhaz[a,2]" %in% variables(fit2))
  expect_range(waic(fit2)$estimates[3, 1], 620, 670)
}))

test_that("ordinal model with grouped thresholds works correctly",
          suppressWarnings({
  set.seed(1234)
  dat <- data.frame(
    y = sample(1:6, 100, TRUE),
    gr = rep(c("a", "b"), each = 50),
    th = rep(5:6, each = 50),
    x = rnorm(100)
  )

  prior <- prior(normal(0,1), class = "Intercept", group = "b")
  fit <- brm(y | thres(th, gr) ~ x, dat, cumulative(), prior = prior)
  print(summary(fit))
  pred <- predict(fit)
  expect_equal(dim(pred), c(nrow(dat), max(dat$th) + 1))
  expect_range(waic(fit)$estimates[3, 1], 350, 400)
  ce <- conditional_effects(fit, categorical = TRUE)
  expect_ggplot(plot(ce, ask = FALSE)[[1]])

  # test incl_thres = TRUE
  thres_minus_eta <- posterior_linpred(fit, incl_thres = TRUE)
  bprep <- prepare_predictions(fit)
  thres <- bprep$thres$thres
  eta <- posterior_linpred(fit)
  gr_unq <- unique(family(fit)$thres$group)
  gr_vec <- fit$data$gr
  nthres_max <- max(
    by(family(fit)$thres, family(fit)$thres$group, function(x) max(x$thres))
  )
  thres_minus_eta_ch <- lapply(setNames(nm = gr_unq), function(gr) {
    thres_gr_nms <- grepl(paste0("^b_Intercept\\[", gr, ","), colnames(thres))
    thres_gr <- thres[, thres_gr_nms]
    eta_gr <- eta[, gr_vec == gr, drop = FALSE]
    thres_minus_eta_ch_gr <- apply(thres_gr, 2, "-", eta_gr)
    thres_minus_eta_ch_gr <- array(
      thres_minus_eta_ch_gr,
      dim = c(nrow(thres_gr), ncol(eta_gr), ncol(thres_gr))
    )
    if (ncol(thres_gr) < nthres_max) {
      dim_NA <- c(
        dim(thres_minus_eta_ch_gr)[-3],
        nthres_max - dim(thres_minus_eta_ch_gr)[3]
      )
      thres_minus_eta_ch_gr <-
        abind::abind(thres_minus_eta_ch_gr, array(dim = dim_NA))
    }
    dimnames(thres_minus_eta_ch_gr) <-
      list(NULL, NULL, as.character(seq_len(nthres_max)))
    return(thres_minus_eta_ch_gr)
  })
  new_arrnms <- dimnames(thres_minus_eta_ch[[1]])
  thres_minus_eta_ch <- abind::abind(thres_minus_eta_ch, along = 2)
  dimnames(thres_minus_eta_ch) <- new_arrnms
  expect_equivalent(thres_minus_eta, thres_minus_eta_ch)
}))

test_that("projpred methods can be run", suppressWarnings({
  fit <- brm(count ~ zAge + zBase + Trt,
             data = epilepsy, family = poisson())
  summary(fit)

  library(projpred)

  # perform variable selection without cross-validation
  vs <- varsel(fit)
  expect_is(vs, "vsel")

  # perform variable selection with cross-validation
  # takes very long and hence commented out
  # cv_vs <- cv_varsel(fit)
  # expect_is(vs, "vsel")
}))

test_that("projpred K-fold reference models use the 'dis' of their refit",
          suppressWarnings({
  fit <- brm(count ~ zAge + zBase, data = epilepsy, refresh = 0)
  # the refit without the first 20 rows stands in for a K-fold refit
  fit_fold <- update(fit, newdata = epilepsy[-(1:20), ], refresh = 0)
  sigma_fold <- as.data.frame(fit_fold, variable = "sigma")$sigma

  ref <- projpred::get_refmodel(fit)
  expect_equal(ref$cvrefbuilder(fit_fold)$dis, sigma_fold)

  dis_user <- rep(1, ndraws(fit))
  ref_user <- projpred::get_refmodel(fit, dis = dis_user)
  expect_equal(ref_user$cvrefbuilder(fit_fold)$dis, dis_user)
}))

test_that("projpred K-fold CV uses a user-supplied 'cvrefbuilder'",
          suppressWarnings({
  fit <- brm(count ~ zAge + zBase, data = epilepsy, refresh = 0)
  # the refit without the first 20 rows stands in for a K-fold refit
  fit_fold <- update(fit, newdata = epilepsy[-(1:20), ], refresh = 0)
  sigma_fold <- as.data.frame(fit_fold, variable = "sigma")$sigma

  cvrb_ks <- integer()
  cvrb <- function(cvfit) {
    cvrb_ks <<- c(cvrb_ks, cvfit$projpred_k)
    sigma <- as.data.frame(cvfit, variable = "sigma")$sigma
    projpred::get_refmodel(cvfit, dis = 2 * sigma,
                           called_from_cvrefbuilder = TRUE)
  }
  ref <- projpred::get_refmodel(fit, cvrefbuilder = cvrb)
  expect_identical(ref$cvrefbuilder, cvrb)
  expect_equal(ref$cvrefbuilder(fit_fold)$dis, 2 * sigma_fold)

  # the 'cvfun' returns the full-data fit for each fold, so no refit is needed
  cvfun <- function(folds, ...) rep(list(fit), max(folds))
  ref <- projpred::get_refmodel(fit, cvfun = cvfun, cvrefbuilder = cvrb)
  cvrb_ks <- integer()
  cvvs <- projpred::cv_varsel(ref, cv_method = "kfold", K = 2,
                              nterms_max = 1, verbose = FALSE)
  expect_is(cvvs, "vsel")
  expect_equal(cvrb_ks, 1:2)

  expect_error(projpred::get_refmodel(fit, cvrefbuilder = "a"),
               "'cvrefbuilder' should be a function")
}))

# simulate the data for the tests of the projpred latent projection
sim_latent_data <- function(N = 100) {
  dat <- data.frame(x1 = rnorm(N), x2 = rnorm(N), x3 = rnorm(N))
  eta <- 0.5 * dat$x1 - 0.3 * dat$x2
  ctime <- runif(N, 0.5, 3)
  y_weib <- rweibull(N, shape = 1.5, scale = exp(eta) / gamma(1 + 1 / 1.5))
  dat$y_weib <- pmin(y_weib, ctime)
  dat$c_weib <- as.numeric(y_weib > ctime)
  y_lnorm <- rlnorm(N, meanlog = eta, sdlog = 0.5)
  dat$y_lnorm <- pmin(y_lnorm, ctime)
  dat$c_lnorm <- ifelse(y_lnorm > ctime, "right", "none")
  dat$y_nb <- rnbinom(N, size = 2, mu = exp(1 + eta))
  y_pois <- rpois(N, exp(1 + eta))
  dat$y_pois <- pmin(y_pois, 5)
  dat$c_pois <- as.numeric(y_pois > 5)
  dat$y_bern <- rbinom(N, 1, plogis(eta))
  dat$n_bb <- 20
  dat$y_bb <- rbinom(N, dat$n_bb, rbeta(N, 5 * plogis(eta), 5 * plogis(-eta)))
  dat$y_zi <- ifelse(runif(N) < 0.2, 0, dat$y_bb)
  dat$eta <- eta
  dat
}

# compare the response-scale functions of brms for the latent projection with
# the brms results, and run varsel(), project(), and proj_predict()
# @param cens censoring indicators of the observations or NULL
# @param trials number of trials of all observations or NULL
# @return a list with the reference model 'ref' and the draws 'ppd' of
#   proj_predict()
check_latent_oscale <- function(fit, cens = NULL, trials = NULL, info = NULL) {
  ref <- projpred::get_refmodel(fit, latent = TRUE)
  # with all draws, the functions give the brms results
  ilpreds <- posterior_linpred(fit, transform = TRUE)
  if (!is.null(trials)) {
    # for families with trials, 'latent_ilink' gives 'trials * mu'
    ilpreds <- trials * ilpreds
  }
  expect_equivalent(ref$family$latent_ilink(t(ref$mu), cl_ref = 1),
                    ilpreds, info = info)
  ll <- ref$family$latent_ll_oscale(
    ilpreds, dis = ref$dis, y_oscale = ref$y_oscale, cens = cens,
    cl_ref = seq_len(ndraws(fit))
  )
  expect_equivalent(ll, log_lik(fit), info = info)
  # projpred versions with argument 'latent_epred' get the expected values
  if ("latent_epred" %in% names(formals(projpred::extend_family))) {
    epred <- ref$family$latent_epred(
      ilpreds, dis = ref$dis, cl_ref = seq_len(ndraws(fit))
    )
    expect_equivalent(epred, posterior_epred(fit), info = info)
  }

  vs <- projpred::varsel(ref, nterms_max = 2, verbose = FALSE)
  smry <- summary(vs, stats = "elpd")$perf_sub
  expect_true(all(is.finite(smry$elpd)), info = info)
  prj <- projpred::project(ref, predictor_terms = "x1", nclusters = 5)
  ppd <- projpred::proj_predict(prj, nresample_clusters = 30)
  expect_equal(dim(ppd), c(30, nobs(fit)), info = info)
  expect_true(all(is.finite(ppd) & ppd >= 0), info = info)
  list(ref = ref, ppd = ppd)
}

test_that("projpred latent projection has response-scale functions",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()

  fits <- list(
    bernoulli = brm(y_bern ~ x1 + x2 + x3, data = dat,
                    family = bernoulli(), refresh = 0),
    poisson = brm(y_pois | cens(c_pois) ~ x1 + x2 + x3, data = dat,
                  family = poisson(), refresh = 0),
    negbinomial = brm(y_nb ~ x1 + x2 + x3, data = dat,
                      family = negbinomial(), refresh = 0),
    weibull = brm(y_weib | cens(c_weib) ~ x1 + x2 + x3, data = dat,
                  family = weibull(), refresh = 0),
    # the censoring times are at most 3, and P(Y > 3) changes the results
    lognormal = brm(y_lnorm | cens(c_lnorm) + trunc(ub = 3) ~ x1 + x2 + x3,
                    data = dat, family = lognormal(), refresh = 0)
  )
  cens <- list(bernoulli = NULL, poisson = dat$c_pois, negbinomial = NULL,
               weibull = dat$c_weib, lognormal = dat$c_lnorm)
  for (fam in names(fits)) {
    fit <- fits[[fam]]
    out <- check_latent_oscale(fit, cens = cens[[fam]], info = fam)
    bounds_attr <- attr(out$ref$family$latent_ll_oscale, "trunc_bounds")
    if (fam == "lognormal") {
      expect_equal(out$ref$dis, as.data.frame(fit, variable = "sigma")$sigma)
      expect_equal(bounds_attr, list(lb = -Inf, ub = 3))
      expect_true(all(out$ppd <= 3))
    } else {
      expect_null(bounds_attr, info = fam)
    }
  }
}))

test_that("projpred latent projection supports families with trials",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()
  N <- nrow(dat)
  trials <- dat$n_bb[1]

  fits <- list(
    beta_binomial = brm(y_bb | trials(n_bb) ~ x1 + x2 + x3, data = dat,
                        family = beta_binomial(), refresh = 0),
    zero_inflated_binomial = brm(y_zi | trials(n_bb) ~ x1 + x2 + x3,
                                 data = dat, refresh = 0,
                                 family = zero_inflated_binomial()),
    zero_inflated_beta_binomial = brm(y_zi | trials(n_bb) ~ x1 + x2 + x3,
                                      data = dat, refresh = 0,
                                      family = zero_inflated_beta_binomial())
  )
  for (fam in names(fits)) {
    fit <- fits[[fam]]
    out <- check_latent_oscale(fit, trials = trials, info = fam)
    ref <- out$ref
    expect_null(attr(ref$family$latent_ll_oscale, "trunc_bounds"), info = fam)
    # the trials go to the response-scale functions, not to the weights
    expect_equal(ref$wobs, rep(1, N), info = fam)
    expect_true(all(out$ppd <= trials & out$ppd == round(out$ppd)),
                info = fam)
    dat_new <- dat[1:5, ]
    expect_equal(ref$extract_model_data(fit, newdata = dat_new,
                                        extract_y = FALSE)$weights,
                 rep(1, 5), info = fam)
    dat_new$n_bb <- 10
    expect_error(
      ref$extract_model_data(fit, newdata = dat_new, extract_y = FALSE),
      paste0("number of trials of the reference model (", trials, ")"),
      fixed = TRUE
    )
  }
}))

test_that("projpred latent functions use refits, clusters, and weights",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()
  N <- nrow(dat)
  fit <- brm(y_weib | cens(c_weib) ~ x1 + x2 + x3, data = dat,
             family = weibull(), refresh = 0)

  # the K-fold reference models use the draws of their refit
  fit_fold <- update(fit, newdata = dat[-(1:20), ], refresh = 0)
  ref_fold <- projpred::get_refmodel(fit, latent = TRUE)$cvrefbuilder(fit_fold)
  ll_fold <- ref_fold$family$latent_ll_oscale(
    posterior_linpred(fit_fold, transform = TRUE), dis = ref_fold$dis,
    y_oscale = ref_fold$y_oscale, cens = dat$c_weib[-(1:20)],
    cl_ref = seq_len(ndraws(fit_fold))
  )
  expect_equivalent(ll_fold, log_lik(fit_fold))

  # clustered draws and observation weights, compared with base R
  ref <- projpred::get_refmodel(fit, latent = TRUE)
  S <- ndraws(fit)
  cl_ref <- rep(1:2, length.out = S)
  wdraws_ref <- seq_len(S)
  shape <- as.data.frame(fit, variable = "shape")$shape
  shape_cl <- c(weighted.mean(shape[cl_ref == 1], wdraws_ref[cl_ref == 1]),
                weighted.mean(shape[cl_ref == 2], wdraws_ref[cl_ref == 2]))
  ilpreds <- rbind(exp(dat$eta), exp(0.5 * dat$eta))
  wobs <- runif(N, 0.5, 2)
  ll <- ref$family$latent_ll_oscale(
    ilpreds, dis = c(NA, NA), y_oscale = dat$y_weib, wobs = wobs,
    cens = dat$c_weib, cl_ref = cl_ref, wdraws_ref = wdraws_ref
  )
  scale <- ilpreds / gamma(1 + 1 / shape_cl)
  y_mat <- matrix(dat$y_weib, 2, N, byrow = TRUE)
  cens_mat <- matrix(dat$c_weib == 1, 2, N, byrow = TRUE)
  ll_base <- ifelse(
    cens_mat,
    pweibull(y_mat, shape_cl, scale, lower.tail = FALSE, log.p = TRUE),
    dweibull(y_mat, shape_cl, scale, log = TRUE)
  )
  expect_equivalent(ll, ll_base * rep(wobs, each = 2))

  idxs <- c(2, 1, 2)
  set.seed(5)
  ppd <- ref$family$latent_ppd_oscale(
    ilpreds[idxs, ], dis_resamp = rep(NA, 3), cl_ref = cl_ref,
    wdraws_ref = wdraws_ref, idxs_prjdraws = idxs
  )
  set.seed(5)
  ppd_base <- rweibull(3 * N, shape_cl[idxs], scale[idxs, ])
  expect_equivalent(ppd, matrix(ppd_base, 3))
}))

test_that("projpred latent projection warns or stops for unsupported models",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()
  N <- nrow(dat)

  # unsupported models give a warning and no functions
  dat$c_int <- ifelse(dat$c_weib == 1, "interval", "none")
  dat$y2 <- dat$y_weib + 1
  dat$ub_var <- 100 + seq_len(N)
  dat$denom <- runif(N, 0.5, 2)
  fits_unsupp <- list(
    "family 'gamma'" = brm(y_weib ~ x1, data = dat, family = Gamma("log"),
                           empty = TRUE),
    "interval censoring" = brm(y_weib | cens(c_int, y2) ~ x1, data = dat,
                               family = weibull(), empty = TRUE),
    "truncation bounds" = brm(y_nb | trunc(ub = ub_var) ~ x1, data = dat,
                              family = negbinomial(), empty = TRUE),
    "addition argument 'rate'" = brm(y_nb | rate(denom) ~ x1, data = dat,
                                     family = negbinomial(), empty = TRUE),
    # brms has no expected values of truncated beta-binomial models
    "addition argument 'trunc' in family 'beta_binomial'" =
      brm(y_bb | trials(n_bb) + trunc(ub = 20) ~ x1, data = dat,
          family = beta_binomial(), empty = TRUE)
  )
  for (unsupp in names(fits_unsupp)) {
    fit_unsupp <- fits_unsupp[[unsupp]]
    expect_warning(
      funs <- brms:::latent_oscale_funs(
        fit_unsupp, family = family(fit_unsupp),
        family_name = family(fit_unsupp)$family,
        bterms = brmsterms(formula(fit_unsupp)), resp = NULL
      ),
      unsupp, fixed = TRUE
    )
    expect_null(funs)
  }

  # the reasons that the warning names
  dat$sdy <- 0.3
  unsupp_reason <- function(formula, family_name) {
    brms:::latent_oscale_unsupported(family_name, brmsterms(formula),
                                     data = dat)
  }
  expect_equal(unsupp_reason(y_weib | mi(sdy) ~ x1, "gaussian"),
               "addition argument 'mi'")
  expect_equal(unsupp_reason(y_weib | trunc(lb = 0) ~ x1, "student"),
               "family 'student' with addition argument(s) 'trunc'")
  expect_null(
    unsupp_reason(y_weib | cens(c_weib) + trunc(lb = 0) ~ x1, "gaussian")
  )

  # weights and trials other than 1 are not supported
  dat$w <- runif(N, 0.5, 2)
  fit_w <- brm(y_weib | weights(w) ~ x1, data = dat, family = weibull(),
               empty = TRUE)
  expect_error(projpred::get_refmodel(fit_w, latent = TRUE),
               "does not support 'weights' or 'trials'")
  fit_tr <- brm(y_nb | trials(200) ~ x1, data = dat, family = binomial(),
                empty = TRUE)
  expect_error(projpred::get_refmodel(fit_tr, latent = TRUE),
               "does not support 'weights' or 'trials'")
}))

test_that("projpred latent projection handles varying numbers of trials",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()
  N <- nrow(dat)
  dat$n_var <- rep(c(10, 20), length.out = N)
  dat$y_var <- pmin(dat$y_bb, dat$n_var)
  fit <- brm(y_var | trials(n_var) ~ x1 + x2, data = dat,
             family = beta_binomial(), refresh = 0)

  # the projection works, and the response-scale results are NA
  expect_warning(
    ref <- projpred::get_refmodel(fit, latent = TRUE),
    "numbers of trials that vary between observations", fixed = TRUE
  )
  expect_equal(ref$wobs, rep(1, N))
  ll <- ref$family$latent_ll_oscale(
    t(ref$mu), dis = ref$dis, y_oscale = ref$y_oscale, cens = NULL,
    cl_ref = seq_len(ndraws(fit))
  )
  expect_true(all(is.na(ll)))
  vs <- projpred::varsel(ref, nterms_max = 1, verbose = FALSE)
  expect_is(vs, "vsel")
  # new data can have other numbers of trials
  dat_new <- dat[1:4, ]
  dat_new$n_var <- 30
  expect_equal(ref$extract_model_data(fit, newdata = dat_new,
                                      extract_y = FALSE)$weights,
               rep(1, 4))

  # a K-fold refit with one number of trials uses the decision of the full
  # fit, so it gets no brms functions and accepts other numbers of trials
  fit_fold <- update(fit, newdata = dat[dat$n_var == 20, ], refresh = 0)
  ref_fold <- ref$cvrefbuilder(fit_fold)
  expect_equal(ref_fold$extract_model_data(fit_fold, newdata = dat[1:4, ],
                                           extract_y = FALSE)$weights,
               rep(1, 4))
  ll_fold <- ref_fold$family$latent_ll_oscale(
    t(ref_fold$mu), dis = ref_fold$dis, y_oscale = ref_fold$y_oscale,
    cens = NULL, cl_ref = seq_len(ndraws(fit_fold))
  )
  expect_true(all(is.na(ll_fold)))
}))

test_that("projpred latent projection keeps user-supplied functions",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()
  N <- nrow(dat)
  fit <- brm(y_weib | cens(c_weib) ~ x1 + x2 + x3, data = dat,
             family = weibull(), refresh = 0)

  # functions passed by the user take precedence
  ll_user <- function(ilpreds, dis, y_oscale, wobs = rep(1, ncol(ilpreds)),
                      cens, cl_ref, wdraws_ref = rep(1, length(cl_ref))) {
    matrix(0, nrow(ilpreds), ncol(ilpreds))
  }
  ref <- projpred::get_refmodel(fit, latent = TRUE, latent_ll_oscale = ll_user)
  expect_identical(ref$family$latent_ll_oscale, ll_user)
  # brms adds none of its functions if the user passes one of them
  ilink_user <- function(lpreds, cl_ref, wdraws_ref = rep(1, length(cl_ref))) {
    exp(lpreds)
  }
  ref <- projpred::get_refmodel(fit, latent = TRUE, latent_ilink = ilink_user)
  S <- ndraws(fit)
  ilpreds <- rbind(exp(dat$eta), exp(0.5 * dat$eta))
  ll <- ref$family$latent_ll_oscale(
    ilpreds, dis = c(NA, NA), y_oscale = dat$y_weib, cens = dat$c_weib,
    cl_ref = rep(1:2, length.out = S), wdraws_ref = seq_len(S)
  )
  expect_true(all(is.na(ll)))

  # with 'latent_y_unqs', projpred uses its functions for categories
  dat$y_fac <- factor(dat$y_bern, labels = c("a", "b"))
  fit_fac <- brm(y_fac ~ x1 + x2 + x3, data = dat, family = bernoulli(),
                 refresh = 0)
  ref <- projpred::get_refmodel(fit_fac, latent = TRUE,
                                latent_y_unqs = c("a", "b"))
  ilpreds <- ref$family$latent_ilink(t(ref$mu),
                                     cl_ref = seq_len(ndraws(fit_fac)))
  expect_equal(dim(ilpreds), c(ndraws(fit_fac), N, 2))
}))

test_that("projpred latent functions work with 'resp' and other links",
          suppressWarnings({
  set.seed(1234)
  dat <- sim_latent_data()
  N <- nrow(dat)
  # projpred versions with argument 'latent_epred' get the expected values
  has_epred <- "latent_epred" %in% names(formals(projpred::extend_family))

  # in multivariate models, the functions are defined for the response 'resp'
  dat$y_g <- rnorm(N, dat$eta)
  fit_mv <- brm(
    bf(y_weib | cens(c_weib) ~ x1 + x2 + x3, family = weibull()) +
      bf(y_nb ~ x1 + x2 + x3, family = negbinomial()) +
      bf(y_g ~ x1 + x2 + x3, family = gaussian()) + set_rescor(FALSE),
    data = dat, refresh = 0
  )
  cens <- list(yweib = dat$c_weib, ynb = NULL, yg = NULL)
  for (resp in names(cens)) {
    ref <- projpred::get_refmodel(fit_mv, resp = resp, latent = TRUE)
    if (resp == "yg") {
      sigma_yg <- as.data.frame(fit_mv, variable = "sigma_yg")$sigma_yg
      expect_equal(ref$dis, sigma_yg)
    }
    ilpreds <- posterior_linpred(fit_mv, resp = resp, transform = TRUE)
    ll <- ref$family$latent_ll_oscale(
      ilpreds, dis = ref$dis, y_oscale = ref$y_oscale, cens = cens[[resp]],
      cl_ref = seq_len(ndraws(fit_mv))
    )
    expect_equivalent(ll, log_lik(fit_mv, resp = resp), info = resp)
    if (has_epred) {
      epred <- ref$family$latent_epred(
        ilpreds, dis = ref$dis, cl_ref = seq_len(ndraws(fit_mv))
      )
      expect_equivalent(epred, posterior_epred(fit_mv, resp = resp),
                        info = resp)
    }
  }

  # with other links than the identity, 'sigma' is not the SD of the latent
  # Gaussian distribution and comes from the draws of the model
  dat$y_lninv <- rlnorm(N, meanlog = 1 / (2 + 0.3 * dat$x1), sdlog = 0.3)
  fit_inv <- brm(y_lninv ~ x1, data = dat,
                 family = lognormal(link = "inverse"), refresh = 0)
  ref <- projpred::get_refmodel(fit_inv, latent = TRUE)
  expect_true(all(is.na(ref$dis)))
  ll <- ref$family$latent_ll_oscale(
    posterior_linpred(fit_inv, transform = TRUE), dis = ref$dis,
    y_oscale = ref$y_oscale, cens = NULL, cl_ref = seq_len(ndraws(fit_inv))
  )
  expect_equivalent(ll, log_lik(fit_inv))
  if (has_epred) {
    epred <- ref$family$latent_epred(
      posterior_linpred(fit_inv, transform = TRUE), dis = ref$dis,
      cl_ref = seq_len(ndraws(fit_inv))
    )
    expect_equivalent(epred, posterior_epred(fit_inv))
  }
}))

test_that("projpred warns about addition terms if 'latent = FALSE'",
          suppressWarnings({
  set.seed(1234)
  N <- 100
  dat <- data.frame(x1 = rnorm(N))
  dat$y <- abs(rnorm(N, dat$x1)) + 0.1
  dat$se_y <- runif(N, 0.5, 1)
  dat$c_y <- ifelse(dat$y > 1, "right", "none")
  dat$c_none <- "none"
  dat$lb_var <- -seq_len(N)
  dat$count <- rpois(N, 3)
  dat$expo <- runif(N, 1, 2)
  dat$sdy <- 0.3

  # the messages of the warnings for a model
  ad_warning <- function(fit) {
    capture_warnings(brms:::warn_projpred_ignored_ad(
      family(fit)$family, brmsterms(formula(fit)), data = fit$data
    ))
  }
  # the terms that the warning names, and if it recommends 'latent = TRUE'
  cases <- list(
    list(fit = brm(y | se(se_y, sigma = TRUE) ~ x1, dat, empty = TRUE),
         terms = "'se'", latent_ok = FALSE),
    list(fit = brm(y | trunc(lb = 0) ~ x1, dat, empty = TRUE),
         terms = "'trunc'", latent_ok = TRUE),
    list(fit = brm(y | trunc(lb = lb_var) ~ x1, dat, empty = TRUE),
         terms = "'trunc'", latent_ok = FALSE),
    list(fit = brm(y | cens(c_y) ~ x1, dat, empty = TRUE),
         terms = "'cens'", latent_ok = TRUE),
    list(fit = brm(y | cens(c_y) + trunc(lb = 0) ~ x1, dat, empty = TRUE),
         terms = "'trunc', 'cens'", latent_ok = TRUE),
    list(fit = brm(count | rate(expo) ~ x1, dat, family = poisson(),
                   empty = TRUE),
         terms = "'rate'", latent_ok = FALSE),
    list(fit = brm(y | mi(sdy) ~ x1, dat, empty = TRUE),
         terms = "'mi'", latent_ok = FALSE)
  )
  for (case in cases) {
    msg <- ad_warning(case$fit)
    expect_true(grepl(paste0("addition argument(s) ", case$terms, " of"), msg,
                      fixed = TRUE), info = case$terms)
    expect_equal(grepl("Use 'latent = TRUE'", msg, fixed = TRUE),
                 case$latent_ok, info = case$terms)
  }
  # no warning without censored observations or addition terms
  fit_none <- brm(y | cens(c_none) ~ x1, dat, empty = TRUE)
  expect_equal(ad_warning(fit_none), character(0))
  expect_equal(ad_warning(brm(y ~ x1, dat, empty = TRUE)), character(0))

  # no warning before the projpred error for an unsupported family
  fit_st <- brm(y | trunc(lb = 0) ~ x1, dat, family = student(), empty = TRUE)
  w <- capture_warnings(
    expect_error(projpred::get_refmodel(fit_st), "not supported")
  )
  expect_false(any(grepl("Projpred ignores", w, fixed = TRUE)))
}))

test_that("projpred stops for 'se' with 'sigma = FALSE' and missing responses",
          suppressWarnings({
  set.seed(1234)
  N <- 100
  dat <- data.frame(x1 = rnorm(N))
  dat$y <- abs(rnorm(N, dat$x1)) + 0.1
  dat$se_y <- runif(N, 0.5, 1)
  dat$y_mis <- dat$y
  dat$y_mis[1:3] <- NA

  fit_se0 <- brm(y | se(se_y) ~ x1, dat, empty = TRUE)
  expect_error(projpred::get_refmodel(fit_se0),
               "'se' with 'sigma = FALSE'", fixed = TRUE)
  fit_mi <- brm(y_mis | mi() ~ x1, dat, empty = TRUE)
  expect_error(projpred::get_refmodel(fit_mi), "missing response values",
               fixed = TRUE)
  expect_error(projpred::get_refmodel(fit_mi, latent = TRUE),
               "missing response values", fixed = TRUE)
}))

test_that("projpred latent projection supports Gaussian 'cens' and 'trunc'",
          suppressWarnings({
  set.seed(1234)
  N <- 100
  dat <- data.frame(x1 = rnorm(N), x2 = rnorm(N), x3 = rnorm(N))
  eta <- 0.5 + 0.5 * dat$x1 - 0.3 * dat$x2
  # normal distribution truncated at 0 and censored at 2
  y <- qnorm(runif(N, pnorm(0, eta), 1), eta)
  dat$c_y <- ifelse(y > 2, "right", "none")
  dat$y <- pmin(y, 2)
  fit <- brm(y | cens(c_y) + trunc(lb = 0) ~ x1 + x2 + x3, data = dat,
             refresh = 0)

  # 'latent = FALSE' gives one warning, also with K-fold CV
  expect_warning(
    projpred::get_refmodel(fit),
    "addition argument(s) 'trunc', 'cens' of the response", fixed = TRUE
  )
  cvfun <- function(folds, ...) rep(list(fit), max(folds))
  w <- capture_warnings(
    projpred::cv_varsel(fit, cv_method = "kfold", K = 2, cvfun = cvfun,
                        nterms_max = 1, verbose = FALSE)
  )
  expect_equal(sum(grepl("Projpred ignores", w, fixed = TRUE)), 1)

  # with 'latent = TRUE' and all draws, the functions give the brms results
  has_epred <- "latent_epred" %in% names(formals(projpred::extend_family))
  ref <- projpred::get_refmodel(fit, latent = TRUE)
  expect_equal(ref$dis, as.data.frame(fit, variable = "sigma")$sigma)
  ilpreds <- posterior_linpred(fit, transform = TRUE)
  ll <- ref$family$latent_ll_oscale(
    ilpreds, dis = ref$dis, y_oscale = ref$y_oscale, cens = dat$c_y,
    cl_ref = seq_len(ndraws(fit))
  )
  expect_equivalent(ll, log_lik(fit))
  if (has_epred) {
    epred <- ref$family$latent_epred(
      ilpreds, dis = ref$dis, cl_ref = seq_len(ndraws(fit))
    )
    expect_equivalent(epred, posterior_epred(fit))
  }
  vs <- projpred::varsel(ref, nterms_max = 2, verbose = FALSE)
  for (oscale in c(TRUE, FALSE)) {
    smry <- summary(vs, stats = "elpd", resp_oscale = oscale)$perf_sub
    expect_true(all(is.finite(smry$elpd)), info = paste("oscale:", oscale))
  }
  prj <- projpred::project(ref, predictor_terms = "x1", nclusters = 5)
  ppd <- projpred::proj_predict(prj, nresample_clusters = 30)
  expect_true(all(is.finite(ppd) & ppd >= 0))
}))

test_that("projpred latent projection warns for 'se'", suppressWarnings({
  set.seed(1234)
  N <- 100
  dat <- data.frame(x1 = rnorm(N))
  dat$se_y <- runif(N, 0.2, 0.5)
  dat$y <- rnorm(N, dat$x1, sqrt(1 + dat$se_y^2))

  fit_se <- brm(y | se(se_y, sigma = TRUE) ~ x1, data = dat, refresh = 0)
  expect_warning(
    ref <- projpred::get_refmodel(fit_se, latent = TRUE),
    "addition argument 'se'", fixed = TRUE
  )
  # with 'se(sigma = TRUE)', 'sigma' is the SD of the latent Gaussian
  # distribution
  expect_equal(ref$dis, as.data.frame(fit_se, variable = "sigma")$sigma)
  # no warning if the user passes the functions
  ll_user <- function(ilpreds, dis, y_oscale, wobs = rep(1, ncol(ilpreds)),
                      cens, cl_ref, wdraws_ref = rep(1, length(cl_ref))) {
    matrix(0, nrow(ilpreds), ncol(ilpreds))
  }
  w <- capture_warnings(
    projpred::get_refmodel(fit_se, latent = TRUE, latent_ll_oscale = ll_user)
  )
  expect_false(any(grepl("not yet implemented in brms", w, fixed = TRUE)))
}))

test_that("projpred latent projection warns for binomial 'trunc'",
          suppressWarnings({
  set.seed(1234)
  N <- 100
  dat <- data.frame(x1 = rnorm(N))
  dat$y_bin <- rbinom(N, 1, plogis(dat$x1))
  # the projpred functions for binomial models ignore truncation
  fit_bin <- brm(y_bin | trials(1) + trunc(ub = 1) ~ x1, data = dat,
                 family = binomial(), refresh = 0)
  expect_warning(
    projpred::get_refmodel(fit_bin, latent = TRUE),
    "'trunc' in family 'binomial'", fixed = TRUE
  )
}))

test_that("projpred methods can be run with monotonic terms",
          suppressWarnings({
  set.seed(1234)
  dat <- data.frame(
    x = sample(0:4, 100, TRUE), z = rnorm(100),
    f = factor(sample(c("lo", "mid", "hi"), 100, TRUE),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
  dat$y <- c(0, 1, 1.5, 1.7, 1.8)[dat$x + 1] +
    c(0, 0.4, 0.6)[as.integer(dat$f)] + 0.5 * dat$z + rnorm(100, sd = 0.5)
  fit <- brm(y ~ mo(x) + mo(f) + z, data = dat, refresh = 0)

  vs <- projpred::varsel(fit)
  expect_is(vs, "vsel")
  expect_setequal(projpred::ranking(vs)$fulldata, c("mo(x)", "mo(f)", "z"))
  prj <- projpred::project(fit, predictor_terms = c("mo(x)", "mo(f)", "z"),
                           ndraws = 25)
  pred_new <- projpred::proj_linpred(prj, newdata = dat[1:5, ])$pred
  pred_train <- projpred::proj_linpred(prj)$pred
  expect_equal(unname(pred_new), unname(pred_train[, 1:5]))

  # me() and mi() terms are not supported
  fit_me <- brm(y ~ me(z, 0.1) + mo(x), data = dat, empty = TRUE)
  expect_error(projpred::get_refmodel(fit_me), "me\\(\\) or mi\\(\\)")

  # the checks also apply to categorical models; 3 categories give the dpars
  # 'mu2' and 'mu3' and no dpar 'mu'
  dat$ycat <- factor(sample(1:3, 100, TRUE))
  fit_me_cat <- brm(ycat ~ me(z, 0.1) + mo(x), data = dat,
                    family = categorical(), empty = TRUE)
  expect_error(projpred::get_refmodel(fit_me_cat), "me\\(\\) or mi\\(\\)")
  fit_gp_cat <- brm(ycat ~ gp(z), data = dat,
                    family = categorical(), empty = TRUE)
  expect_error(projpred::get_refmodel(fit_gp_cat),
               "only supports standard multilevel")
}))

test_that("alternative algorithms can be used", suppressWarnings({
  fit <- brm(
    count ~ zBase, data = epilepsy,
    backend = "cmdstanr", algorithm = "meanfield"
  )
  summary(fit)
  expect_is(fit, "brmsfit")

  fit <- brm(
    count ~ zBase, data = epilepsy,
    backend = "cmdstanr", algorithm = "pathfinder"
  )
  summary(fit)
  expect_is(fit, "brmsfit")

  fit <- brm(
    count ~ zBase, data = epilepsy,
    backend = "cmdstanr", algorithm = "laplace"
  )
  summary(fit)
  expect_is(fit, "brmsfit")
}))

test_that(paste(
  "Families sratio() and cratio() are equivalent for symmetric distribution",
  "functions (here only testing the logit link)"
), suppressWarnings({
  set.seed(1234)
  dat2 <- data.frame(
    rating = sample(1:4, 50, TRUE),
    subject = rep(1:10, 5),
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50)
  )
  warmup <- 150
  iter <- 200
  chains <- 1

  fit_sratio <- brm(
    bf(rating ~ x1 + cs(x2) + (cs(x2)||subject), disc ~ 1),
    data = dat2, family = sratio(),
    warmup = warmup, iter = iter, chains = chains,
    seed = 533273
  )
  draws_sratio <- as.matrix(fit_sratio)

  fit_cratio <- brm(
    bf(rating ~ x1 + cs(x2) + (cs(x2)||subject), disc ~ 1),
    data = dat2, family = cratio(),
    warmup = warmup, iter = iter, chains = chains,
    seed = 533273
  )
  draws_cratio <- as.matrix(fit_cratio)

  expect_equal(draws_sratio, draws_cratio)
}))


test_that("prior tags correctly go to priorsense data", {
  library(priorsense)

  warmup <- 150
  iter <- 200
  chains <- 1

  set.seed(1234)
  N <- 40
  dat <- data.frame(
    count = rpois(N, lambda = 20),
    visit = factor(rep(1:4, each = N/4)),
    patient = factor(rep(1:(N/4), 4)),
    Trt = factor(sample(0:1, N, TRUE))
  )

  fit_tags <- SW(brm(
    formula = count ~ Trt + (1 | patient) + (1 + Trt | visit),
    data = dat,
    warmup = warmup, iter = iter, chains = chains,
    prior = c(prior(normal(0, 1), class = sd, tag = "prior_tag1"),
              prior(normal(0, 5), class = b, tag = "prior_tag2"),
              prior(normal(0, 0.5), coef = "Trt1", tag = "prior_tag3"),
              prior(normal(0, 10), class = "Intercept", tag = "prior_tag4"),
              prior(lkj_corr_cholesky(3), class = "L", group = "visit", tag = "prior_tag5")
              ),
    family = poisson()
  ))

  psd <- create_priorsense_data(fit_tags)

  expect_setequal(
    variables(psd$log_prior),
    c(
      "lprior",
      "lprior_prior_tag1",
      "lprior_prior_tag2",
      "lprior_prior_tag3",
      "lprior_prior_tag4",
      "lprior_prior_tag5"
    )
  )
})
