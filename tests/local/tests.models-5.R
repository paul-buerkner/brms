source("setup_tests_local.R")

test_that("global shrinkage priors work correctly", {
  set.seed(2563)
  dat <- epilepsy
  dat$x1 <- sample(1:4, nrow(dat), TRUE)
  dat$x2 <- sample(1:10, nrow(dat), TRUE)
  dat$x3 <- sample(1:10, nrow(dat), TRUE)

  bprior <- prior(R2D2(), class = b) +
    prior(R2D2(main = FALSE), class = sd) +
    prior(R2D2(main = FALSE), class = sds) +
    prior(R2D2(main = FALSE), class = sdgp) +
    prior(R2D2(main = FALSE), class = ar) +
    prior(R2D2(main = FALSE), class = ma)
  bform <- bf(count ~ Trt * Base + Age + mo(x1) + (1|patient) +
                gp(x2) + s(x3) + arma(p = 2, q = 2, gr = patient))
  fit <- brm(bform, data = dat, prior = bprior, cores = 4,
             control = list(adapt_delta = 0.95), seed = 8892)

  classes <- c("sdb", "sdbsp", "sdbs", "sdar", "sdma")
  for (cl in classes) {
    expect_true(any(grepl(paste0("^", cl), variables(fit))))
  }
  expect_range(SW(waic(fit))$estimates[3, 1], 1580, 1660)
})

test_that("Addition argument 'subset' works correctly", {
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
})

test_that("Cox models work correctly", {
  set.seed(12345)
  covs <- data.frame(id  = 1:200, trt = stats::rbinom(200, 1L, 0.5))
  d1 <- simsurv::simsurv(lambdas = 0.1, gammas  = 1.5, betas = c(trt = -0.5),
                         x = covs, maxt  = 5)
  d1 <- merge(d1, covs)

  fit1 <- brm(eventtime | cens(1 - status) ~ 1 + trt,
              data = d1, family = brmsfamily("cox"), refresh = 0)
  print(summary(fit1))
  expect_range(posterior_summary(fit1)["b_trt", "Estimate"], -0.70, -0.30)
  expect_range(waic(fit1)$estimates[3, 1], 620, 670)
})

test_that("ordinal model with grouped thresholds works correctly", {
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
})



test_that("projpred methods can be run", {
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
})

test_that(paste(
  "Families sratio() and cratio() are equivalent for symmetric distribution",
  "functions (here only testing the logit link)"
), {
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
})
