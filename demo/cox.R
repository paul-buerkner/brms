library(brms)

# These are demonstration functions, not part of the brms API. Their use of
# brms::: helpers and prepared-object fields is implementation-dependent.
# Baseline functions return draws by times; conditional functions return
# draws by observations and describe untruncated event times.
cox_bhaz <- function(t, object, stratum = 1L, ...) {
  prep <- brms:::prepare_cox_predictions(object, ...)
  stratum <- brms:::as_one_integer(stratum)
  if (stratum < 1L || stratum > length(prep$bhaz$sbhaz)) {
    stop("Invalid baseline hazard stratum.")
  }
  brms:::bhaz_at(t, prep$bhaz$sbhaz[[stratum]], prep$bhaz$basis,
                 prep$bhaz$tail, paired = FALSE)
}

cox_cumbhaz <- function(t, object, stratum = 1L, ...) {
  prep <- brms:::prepare_cox_predictions(object, ...)
  stratum <- brms:::as_one_integer(stratum)
  if (stratum < 1L || stratum > length(prep$bhaz$sbhaz)) {
    stop("Invalid baseline hazard stratum.")
  }
  brms:::bhaz_at(t, prep$bhaz$sbhaz[[stratum]], prep$bhaz$basis,
                 prep$bhaz$tail, integrate = TRUE, paired = FALSE)
}

cox_quantity <- function(t, prep, integrate, sort) {
  if (!length(t) %in% c(1L, prep$nobs)) {
    stop("'t' must have length 1 or the number of observations.")
  }
  t <- rep_len(t, prep$nobs)
  if (!sort && length(prep$old_order)) t <- t[order(prep$old_order)]
  out <- matrix(NA_real_, prep$ndraws, prep$nobs)
  for (s in seq_along(prep$bhaz$sbhaz)) {
    take <- which(prep$bhaz$Jgrbhaz == s)
    if (!length(take)) next
    out[, take] <- brms:::bhaz_at(
      t[take], prep$bhaz$sbhaz[[s]], prep$bhaz$basis, prep$bhaz$tail,
      integrate = integrate, paired = FALSE
    ) * prep$dpars$mu[, take, drop = FALSE]
  }
  brms:::reorder_obs(out, prep$old_order, sort = sort)
}

cox_hazard <- function(t, object, sort = FALSE, ...) {
  cox_quantity(t, brms:::prepare_cox_predictions(object, ...), FALSE, sort)
}

cox_cumhaz <- function(t, object, sort = FALSE, ...) {
  cox_quantity(t, brms:::prepare_cox_predictions(object, ...), TRUE, sort)
}

cox_surv <- function(t, object, sort = FALSE, ...) {
  exp(-cox_cumhaz(t, object, sort = sort, ...))
}

# Full curves: draws by times by observations.
# Streaming summaries: times by quantiles by observations.
cox_surv_curve <- function(times, object, sort = FALSE, ...) {
  prep <- brms:::prepare_cox_predictions(object, ...)
  cox_curve(times, prep, sort = sort)
}

cox_surv_summary <- function(times, object, probs = c(0.025, 0.5, 0.975),
                              sort = FALSE, ...) {
  if (!length(probs) || anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop("'probs' must be a non-empty vector of probabilities.")
  }
  prep <- brms:::prepare_cox_predictions(object, ...)
  cox_curve(times, prep, probs = probs, sort = sort)
}

cox_curve <- function(times, prep, probs = NULL, sort = FALSE) {
  if (!length(times) || anyNA(times) || any(!is.finite(times) | times < 0)) {
    stop("'times' must be a non-empty vector of finite nonnegative times.")
  }
  if (is.null(probs)) {
    out <- array(NA_real_, c(prep$ndraws, length(times), prep$nobs))
  } else {
    out <- array(NA_real_, c(length(probs), length(times), prep$nobs))
  }
  for (s in seq_along(prep$bhaz$sbhaz)) {
    H0 <- brms:::bhaz_at(times, prep$bhaz$sbhaz[[s]], prep$bhaz$basis,
                        prep$bhaz$tail, integrate = TRUE, paired = FALSE)
    for (i in which(prep$bhaz$Jgrbhaz == s)) {
      S <- exp(-H0 * prep$dpars$mu[, i])
      out[, , i] <- if (is.null(probs)) S else
        t(matrixStats::colQuantiles(S, probs = probs, drop = FALSE))
    }
  }
  # reorder_obs expects observations in the second dimension.
  out <- aperm(brms:::reorder_obs(aperm(out, c(1, 3, 2)), prep$old_order, sort),
               c(1, 3, 2))
  if (!is.null(probs)) out <- aperm(out, c(2, 1, 3))
  out
}

# Censoring is applied explicitly, never inferred by posterior_predict().
# cens_times can be a scalar, one value per observation, a matching matrix,
# or "reverse_km"; Inf means no additional censoring.
cox_censor <- function(yrep, cens_times, y = NULL, status = NULL) {
  if (!is.matrix(yrep) || !is.numeric(yrep) || anyNA(yrep)) {
    stop("'yrep' must be a numeric matrix without missing values.")
  }
  if (is.character(cens_times)) {
    match.arg(cens_times, "reverse_km")
    ct <- matrix(cox_censor_times(y, status, length(yrep)), nrow(yrep))
  } else if (is.matrix(cens_times)) {
    if (!identical(dim(cens_times), dim(yrep))) {
      stop("'cens_times' must have the same dimensions as 'yrep'.")
    }
    ct <- cens_times
  } else {
    if (!length(cens_times) %in% c(1L, ncol(yrep))) {
      stop("'cens_times' must have length 1 or the number of observations.")
    }
    ct <- matrix(cens_times, nrow(yrep), ncol(yrep), byrow = TRUE)
  }
  if (!is.numeric(ct) || anyNA(ct) || any(ct <= 0)) {
    stop("Censoring times must be positive, with Inf for no censoring.")
  }
  time <- pmin(abs(yrep), ct)
  ifelse(yrep < 0 | abs(yrep) > ct, -time, time)
}

# Reverse-KM simulation uses binary event indicators, not brms censoring codes.
cox_censor_times <- function(y, status, n = length(y)) {
  y <- as.numeric(y)
  status <- as.numeric(status)
  n <- brms:::as_one_integer(n)
  if (!length(y) || length(y) != length(status) || anyNA(status) ||
      any(!status %in% 0:1) || any(!is.finite(y) | y <= 0) || n < 0L) {
    stop("Supply positive finite 'y', matching binary 'status', and n >= 0.")
  }
  if (all(status == 1)) return(rep(Inf, n))
  brms:::require_package("prodlim")
  dat <- data.frame(time = y, status = status)
  fit <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = dat,
                          reverse = TRUE)
  jumps <- sort(unique(y))
  G <- as.vector(predict(fit, times = jumps))
  if (anyNA(G)) {
    stop("The reverse Kaplan-Meier estimate is undefined at observed times.")
  }
  index <- findInterval(-runif(n), -G) + 1L
  out <- jumps[pmin(index, length(jumps))]
  out[index > length(jumps)] <- Inf
  out
}

cox_decode <- function(x) {
  list(time = abs(x), status = (x > 0) * 1L)
}

# Running this demo fits a model using the configured brms Stan backend.
set.seed(123)
x <- rnorm(80)
event_time <- rexp(length(x), rate = exp(0.5 * x) / 5)
censor_time <- rexp(length(x), rate = 1 / 10)
dat <- data.frame(
  time = pmin(event_time, censor_time),
  censored = as.integer(event_time > censor_time),
  x = x
)
fit_cox <- brm(
  time | cens(censored) ~ x, data = dat, family = cox(),
  prior = prior(normal(0, 1), class = "b"),
  chains = 2, cores = 1, iter = 1000, seed = 123, refresh = 0
)

# Prepare once; all subsequent calls use the same posterior draws and tail.
# Constant and loglinear tails are extrapolation assumptions, not observations.
newdata <- data.frame(x = c(-1, 0, 1))
prep <- prepare_predictions(
  fit_cox, newdata = newdata, draw_ids = 1:100,
  check_response = FALSE, tail = "constant"
)
times <- seq(0, max(dat$time), length.out = 30)
horizon <- median(dat$time)
baseline_hazard <- cox_bhaz(times, prep)
baseline_cumulative_hazard <- cox_cumbhaz(times, prep)
hazard <- cox_hazard(horizon, prep)
cumulative_hazard <- cox_cumhaz(horizon, prep)
survival <- cox_surv(horizon, prep)
curves <- cox_surv_curve(times, prep)
curve_summary <- cox_surv_summary(times, prep)
print(curve_summary[, , 2])

# Native predictions are event times; RMST caps times at the explicit horizon.
# Predictive CDFs and RMST also condition on model truncation when present,
# unlike the untruncated survival recipes above.
yrep <- posterior_predict(prep)
rmst <- posterior_epred(
  prep, horizon = horizon, dpar = NULL, nlpar = NULL, sort = FALSE
)
print(posterior_summary(rmst))

# Administrative censoring preserves the draws-by-observations shape.
censored_yrep <- cox_censor(yrep, horizon)
decoded <- cox_decode(censored_yrep)
print(colMeans(decoded$status))

# With no tail, negative draws encode censoring at the spline boundary.
# Do not apply ordinary event-time summaries to these signed values.
yrep_none <- posterior_predict(
  fit_cox, newdata = newdata, draw_ids = 1:100, tail = "none"
)
decoded_none <- cox_decode(yrep_none)
print(colMeans(decoded_none$status))

if (requireNamespace("prodlim", quietly = TRUE)) {
  censoring_times <- cox_censor_times(dat$time, 1L - dat$censored, n = 100)
  observed_prep <- prepare_predictions(fit_cox, draw_ids = 1:100)
  censored_reverse_km <- cox_censor(
    posterior_predict(observed_prep), "reverse_km",
    y = dat$time, status = 1L - dat$censored
  )
  print(colMeans(cox_decode(censored_reverse_km)$status))
} else {
  message("Install 'prodlim' to run the optional reverse-KM censoring example.")
}
