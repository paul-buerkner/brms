params <-
list(EVAL = TRUE)

## ---- SETTINGS-knitr, include=FALSE-----------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "jpeg",
  dpi = 100,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)
library(ggplot2)
library(brms)
theme_set(theme_default())

## ---- fake-data-sim, include=FALSE, eval=TRUE-------------------------------------------
set.seed(54647)
# number of observations
N <- 1E4
# number of group levels
G <- round(N / 10)
# number of predictors
P <- 3
# regression coefficients
beta <- rnorm(P)

# sampled covariates, group means and fake data
fake <- matrix(rnorm(N * P), ncol = P)
dimnames(fake) <- list(NULL, paste0("x", 1:P))

# fixed effect part and sampled group membership
fake <- transform(
  as.data.frame(fake),
  theta = fake %*% beta,
  g = sample.int(G, N, replace=TRUE)
)

# add random intercept by group
fake  <- merge(fake, data.frame(g = 1:G, eta = rnorm(G)), by = "g")

# linear predictor
fake  <- transform(fake, mu = theta + eta)

# sample Poisson data
fake  <- transform(fake, y = rpois(N, exp(mu)))

# shuffle order of data rows to ensure even distribution of computational effort
fake <- fake[sample.int(N, N),]

# drop not needed row names
rownames(fake) <- NULL

## ---- model-poisson, include=FALSE------------------------------------------------------
model_poisson <- brm(
  y ~ 1 + x1 + x2 + (1 | g),
  data = fake,
  family = poisson(),
  iter = 500, # short sampling to speedup example
  chains = 2,
  prior = prior(normal(0,1), class = b) +
    prior(constant(1), class = sd, group = g),
  backend = "cmdstanr",
  threads = threading(4)
)

## ---- benchmark, include=FALSE----------------------------------------------------------
# Benchmarks given model with cross-product of tuning parameters CPU
# cores, grainsize and iterations. Models are run with either static
# or non-static scheduler and initial values are set by default to 0 on the
# unconstrained scale. Function returns a data-frame with the
# cross-product of the tuning parameters and as result column the
# respective runtime.
benchmark_threading <- function(model, cores = 1, grainsize = 1, iter = 100,
                                static = FALSE) {

    winfo <- extract_warmup_info(model)
    sims  <- rstan::extract(model$fit)
    init <- list(extract_draw(sims, 1))

    scaling_model <- update(
        model, refresh = 0,
        threads = threading(1, grainsize = grainsize[1], static = static),
        chains = 1, iter = 2, backend = "cmdstanr"
    )

    run_benchmark <- function(cores, size, iter) {
        bench_fit <- update(
            scaling_model, warmup=0, iter = iter,
            chains = 1, seed = 1234, init = init, refresh = 0, save_warmup=TRUE,
            threads = threading(cores, grainsize = size, static = static),
            inv_metric=winfo$inv_metric[[1]],
            step_size=winfo$step_size[[1]],
            adapt_engaged=FALSE
        )
        lf <- sum(subset(nuts_params(bench_fit, inc_warmup=TRUE), Parameter=="n_leapfrog__")$Value)
        elapsed <- sum(colSums(rstan::get_elapsed_time(bench_fit$fit)))

        c(num_leapfrog=lf, runtime=elapsed)
    }

    cases <- expand.grid(cores = cores, grainsize = grainsize, iter = iter)
    res <- with(cases, mapply(run_benchmark, cores, grainsize, iter))
    cbind(cases, as.data.frame(t(res)))
}

benchmark_reference <- function(model, iter=100, init=0) {
    winfo <- extract_warmup_info(model)
    sims  <- rstan::extract(model$fit)
    init <- list(extract_draw(sims, 1))

    ref_model <- update(
        model, refresh = 0,
        threads = NULL,
        chains = 1, iter = 2, backend = "cmdstanr"
    )

    run_benchmark_ref <- function(iter_bench) {
        bench_fit <- update(
            ref_model, warmup=0, iter = iter_bench,
            chains = 1, seed = 1234, init = init, refresh = 0,
            inv_metric=winfo$inv_metric[[1]],
            step_size=winfo$step_size[[1]],
            adapt_engaged=FALSE
        )

        lf <- sum(subset(nuts_params(bench_fit, inc_warmup=TRUE), Parameter=="n_leapfrog__")$Value)
        elapsed <- sum(colSums(rstan::get_elapsed_time(bench_fit$fit)))

        c(num_leapfrog=lf, runtime=elapsed)
    }

    ref <- sapply(iter, run_benchmark_ref)
    ref <- cbind(as.data.frame(t(ref)), iter=iter)
    ref
}

extract_warmup_info <- function(bfit) {
    adapt  <- lapply(rstan::get_adaptation_info(bfit$fit), strsplit, split="\\n")
    step_size  <- lapply(adapt, function(a) as.numeric(strsplit(a[[1]][[1]], " = ")[[1]][2]))
    inv_metric <- lapply(adapt, function(a) as.numeric(strsplit(sub("^# ", "", a[[1]][[3]]), ", ")[[1]]))
    list(step_size=step_size, inv_metric=inv_metric)
}

extract_draw <- function(sims, draw) lapply(sims, abind::asub, idx=draw, dim=1, drop=FALSE)


## ---- eval=FALSE------------------------------------------------------------------------
#  fit_serial <- brm(
#    count ~ zAge + zBase * Trt + (1|patient),
#    data = epilepsy, family = poisson(),
#    chains = 4, cores = 4, backend = "cmdstanr"
#  )

## ---- eval=FALSE------------------------------------------------------------------------
#  fit_parallel <- update(
#    fit_serial, chains = 2, cores = 2,
#    backend = "cmdstanr", threads = threading(2)
#  )

## ---------------------------------------------------------------------------------------
kable(head(fake, 10), digits = 3)

## ---- eval=FALSE------------------------------------------------------------------------
#  model_poisson <- brm(
#    y ~ 1 + x1 + x2 + (1 | g),
#    data = fake,
#    family = poisson(),
#    iter = 500, # short sampling to speedup example
#    chains = 2,
#    prior = prior(normal(0,1), class = b) +
#      prior(constant(1), class = sd, group = g),
#    backend = "cmdstanr",
#    threads = threading(4)
#  )

## ---- chunking-scale, message=FALSE, warning=FALSE, results='hide'----------------------
chunking_bench <- transform(
    data.frame(chunks = 4^(0:3)),
    grainsize = ceiling(N / chunks)
)

iter_test <- c(10, 20, 40)  # very short test runs
scaling_chunking <- benchmark_threading(
  model_poisson,
  cores = 1,
  grainsize = chunking_bench$grainsize,  # test various grainsizes
  iter = iter_test,
  static = TRUE  # with static partitioner
)

# run as reference the model *without* reduce_sum
ref <- benchmark_reference(model_poisson, iter_test)

# for additional data munging please refer to the appendix

## ---- munge-chunking-scaling, include=FALSE---------------------------------------------
scaling_chunking <- merge(scaling_chunking, chunking_bench, by = "grainsize")

single_chunk  <- transform(
    subset(scaling_chunking, chunks == 1),
    num_leapfrog_single = num_leapfrog, num_leapfrog = NULL,
    runtime_single = runtime, runtime = NULL,
    grainsize = NULL, chunks=NULL
)

scaling_chunking <- transform(
    merge(scaling_chunking, single_chunk),
    slowdown = runtime/runtime_single,
    iter = factor(iter),
    runtime_single = NULL
)

ref <- transform(ref, iter=factor(iter))

## ---------------------------------------------------------------------------------------
ggplot(scaling_chunking) +
    aes(chunks, slowdown, colour = iter, shape = iter) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = scaling_chunking$chunks) +
    scale_y_log10(breaks=seq(0.8, 2.5, by=0.1)) +
    ggtitle("Slowdown with increasing number of chunks")

ggplot(scaling_chunking) +
    aes(chunks, 1E3 * runtime/num_leapfrog, colour = iter, shape=iter) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = scaling_chunking$chunks) +
    scale_y_log10(breaks=seq(0.1, 2.0, by=0.1)) +
    geom_hline(data=ref, aes(yintercept=1E3 * runtime/num_leapfrog, colour=iter), linetype=I(2)) +
    ggtitle("Time per leapfrog step vs number of chunks",
            "Dashed line is reference model without reduce_sum") +
    ylab("Time per leapfrog step [ms]")



## ---- speedup-scale, message=FALSE, warning=FALSE, results='hide'-----------------------
num_cpu <- parallel::detectCores(logical = FALSE)
num_cpu_logical <- parallel::detectCores(logical = TRUE)
grainsize_default <- ceiling(N / (2 * num_cpu))
cores <- c(2^seq(0, floor(log2(num_cpu_logical))), num_cpu, num_cpu_logical)
cores <- sort(unique(cores))
grainsize <- c(grainsize_default, grainsize_default/2, grainsize_default/4)
grainsize <- round(grainsize)

iter_scaling <- 20
scaling_cores <- benchmark_threading(
  model_poisson,
  cores = cores,
  grainsize = grainsize,
  iter = iter_scaling,
  static = FALSE
)

single_core  <- transform(
    subset(scaling_cores, cores == 1),
    runtime_single = runtime,
    num_leapfrog=NULL, runtime=NULL, cores = NULL
)

scaling_cores <- transform(
  merge(scaling_cores, single_core),
  speedup = runtime_single/runtime,
  grainsize = factor(grainsize)
)

## ---------------------------------------------------------------------------------------
ggplot(scaling_cores) +
    aes(cores, runtime, shape = grainsize, color = grainsize) +
    geom_vline(xintercept = num_cpu, linetype = 3) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = scaling_cores$cores) +
    scale_y_log10(breaks=seq(0.1, 1.4, by=0.1)) +
    theme(legend.position = c(0.85, 0.8)) +
    geom_hline(data=subset(ref, iter==iter_scaling), aes(yintercept=runtime), linetype=I(2)) +
    ggtitle("Runtime with varying number of cores",
            "Dashed line is reference model without reduce_sum")

ggplot(scaling_cores) +
  aes(cores, speedup, shape = grainsize, color = grainsize) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_vline(xintercept = num_cpu, linetype = 3) +
  geom_line() + geom_point() +
  scale_x_log10(breaks=scaling_cores$cores) +
  scale_y_log10(breaks=scaling_cores$cores) +
  theme(aspect.ratio = 1) +
  coord_fixed(xlim = c(1, num_cpu_logical), ylim = c(1, num_cpu_logical)) +
  ggtitle("Relative speedup vs 1 core")

## ---------------------------------------------------------------------------------------
kable(scaling_cores, digits = 2)

## ---- eval=FALSE------------------------------------------------------------------------
#  set.seed(54647)
#  # number of observations
#  N <- 1E4
#  # number of group levels
#  G <- round(N / 10)
#  # number of predictors
#  P <- 3
#  # regression coefficients
#  beta <- rnorm(P)
#  
#  # sampled covariates, group means and fake data
#  fake <- matrix(rnorm(N * P), ncol = P)
#  dimnames(fake) <- list(NULL, paste0("x", 1:P))
#  
#  # fixed effect part and sampled group membership
#  fake <- transform(
#    as.data.frame(fake),
#    theta = fake %*% beta,
#    g = sample.int(G, N, replace=TRUE)
#  )
#  
#  # add random intercept by group
#  fake  <- merge(fake, data.frame(g = 1:G, eta = rnorm(G)), by = "g")
#  
#  # linear predictor
#  fake  <- transform(fake, mu = theta + eta)
#  
#  # sample Poisson data
#  fake  <- transform(fake, y = rpois(N, exp(mu)))
#  
#  # shuffle order of data rows to ensure even distribution of computational effort
#  fake <- fake[sample.int(N, N),]
#  
#  # drop not needed row names
#  rownames(fake) <- NULL

## ---- eval=FALSE------------------------------------------------------------------------
#  model_poisson <- brm(
#    y ~ 1 + x1 + x2 + (1 | g),
#    data = fake,
#    family = poisson(),
#    iter = 500, # short sampling to speedup example
#    chains = 2,
#    prior = prior(normal(0,1), class = b) +
#      prior(constant(1), class = sd, group = g),
#    backend = "cmdstanr",
#    threads = threading(4)
#  )

## ---- eval=FALSE------------------------------------------------------------------------
#  # Benchmarks given model with cross-product of tuning parameters CPU
#  # cores, grainsize and iterations. Models are run with either static
#  # or non-static scheduler and initial values are set by default to 0 on the
#  # unconstrained scale. Function returns a data-frame with the
#  # cross-product of the tuning parameters and as result column the
#  # respective runtime.
#  benchmark_threading <- function(model, cores = 1, grainsize = 1, iter = 100,
#                                  static = FALSE) {
#  
#      winfo <- extract_warmup_info(model)
#      sims  <- rstan::extract(model$fit)
#      init <- list(extract_draw(sims, 1))
#  
#      scaling_model <- update(
#          model, refresh = 0,
#          threads = threading(1, grainsize = grainsize[1], static = static),
#          chains = 1, iter = 2, backend = "cmdstanr"
#      )
#  
#      run_benchmark <- function(cores, size, iter) {
#          bench_fit <- update(
#              scaling_model, warmup=0, iter = iter,
#              chains = 1, seed = 1234, init = init, refresh = 0, save_warmup=TRUE,
#              threads = threading(cores, grainsize = size, static = static),
#              inv_metric=winfo$inv_metric[[1]],
#              step_size=winfo$step_size[[1]],
#              adapt_engaged=FALSE
#          )
#          lf <- sum(subset(nuts_params(bench_fit, inc_warmup=TRUE), Parameter=="n_leapfrog__")$Value)
#          elapsed <- sum(colSums(rstan::get_elapsed_time(bench_fit$fit)))
#  
#          c(num_leapfrog=lf, runtime=elapsed)
#      }
#  
#      cases <- expand.grid(cores = cores, grainsize = grainsize, iter = iter)
#      res <- with(cases, mapply(run_benchmark, cores, grainsize, iter))
#      cbind(cases, as.data.frame(t(res)))
#  }
#  
#  benchmark_reference <- function(model, iter=100, init=0) {
#      winfo <- extract_warmup_info(model)
#      sims  <- rstan::extract(model$fit)
#      init <- list(extract_draw(sims, 1))
#  
#      ref_model <- update(
#          model, refresh = 0,
#          threads = NULL,
#          chains = 1, iter = 2, backend = "cmdstanr"
#      )
#  
#      run_benchmark_ref <- function(iter_bench) {
#          bench_fit <- update(
#              ref_model, warmup=0, iter = iter_bench,
#              chains = 1, seed = 1234, init = init, refresh = 0,
#              inv_metric=winfo$inv_metric[[1]],
#              step_size=winfo$step_size[[1]],
#              adapt_engaged=FALSE
#          )
#  
#          lf <- sum(subset(nuts_params(bench_fit, inc_warmup=TRUE), Parameter=="n_leapfrog__")$Value)
#          elapsed <- sum(colSums(rstan::get_elapsed_time(bench_fit$fit)))
#  
#          c(num_leapfrog=lf, runtime=elapsed)
#      }
#  
#      ref <- sapply(iter, run_benchmark_ref)
#      ref <- cbind(as.data.frame(t(ref)), iter=iter)
#      ref
#  }
#  
#  extract_warmup_info <- function(bfit) {
#      adapt  <- lapply(rstan::get_adaptation_info(bfit$fit), strsplit, split="\\n")
#      step_size  <- lapply(adapt, function(a) as.numeric(strsplit(a[[1]][[1]], " = ")[[1]][2]))
#      inv_metric <- lapply(adapt, function(a) as.numeric(strsplit(sub("^# ", "", a[[1]][[3]]), ", ")[[1]]))
#      list(step_size=step_size, inv_metric=inv_metric)
#  }
#  
#  extract_draw <- function(sims, draw) lapply(sims, abind::asub, idx=draw, dim=1, drop=FALSE)
#  

## ---- eval=FALSE------------------------------------------------------------------------
#  scaling_chunking <- merge(scaling_chunking, chunking_bench, by = "grainsize")
#  
#  single_chunk  <- transform(
#      subset(scaling_chunking, chunks == 1),
#      num_leapfrog_single = num_leapfrog, num_leapfrog = NULL,
#      runtime_single = runtime, runtime = NULL,
#      grainsize = NULL, chunks=NULL
#  )
#  
#  scaling_chunking <- transform(
#      merge(scaling_chunking, single_chunk),
#      slowdown = runtime/runtime_single,
#      iter = factor(iter),
#      runtime_single = NULL
#  )
#  
#  ref <- transform(ref, iter=factor(iter))

