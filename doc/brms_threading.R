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
  dev = "png",
  dpi = 150,
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
  iter = 100, # short sampling to speedup example
  prior = prior(normal(0,1), class = b) +
    prior(constant(1), class = sd, group = g),
  backend = "cmdstanr",
  threads = threading(2)
)

## ---- benchmark, include=FALSE----------------------------------------------------------
# Benchmarks given model with cross-product of tuning parameters CPU
# cores, grainsize and iterations. Models are run with either static
# or non-static scheduler and inits is set by default to 0 on the
# unconstrained scale. Function returns a data-frame with the
# cross-product of the tuning parameters and as result column the
# respective runtime.
benchmark_threading <- function(model, cores = 1, grainsize = 1, iter = 100, 
                                static = FALSE, inits = 0) {

  scaling_model <- update(
    model, refresh = 0, 
    threads = threading(1, grainsize = grainsize[1], static = static), 
    chains = 1, iter = 2, backend = "cmdstanr"
  )

  run_benchmark <- function(cores, size, iter) {
    unname(
      system.time(
        update(
          scaling_model, iter = iter, 
          chains = 1, seed = 1234, inits = inits, refresh = 0,
          threads = threading(cores, grainsize = size, static = static)
        )
      )["elapsed"]
    )
  }

  cases <- expand.grid(cores = cores, grainsize = grainsize, iter = iter)
  cases$runtime <- with(cases, mapply(run_benchmark, cores, grainsize, iter))
  cases
}

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
#    iter = 100, # short sampling to speedup example
#    prior = prior(normal(0,1), class = b) +
#      prior(constant(1), class = sd, group = g),
#    backend = "cmdstanr",
#    threads = threading(2)
#  )

## ---- chunking-scale, message=FALSE, warning=FALSE, results='hide'----------------------
chunking_bench <- transform(
  data.frame(chunks = 2^(0:4)), 
  grainsize = ceiling(N / chunks)
)

scaling_chunking <- benchmark_threading(
  model_poisson,
  cores = 1,                         
  grainsize = chunking_bench$grainsize,  # test various grainsizes
  iter = c(25, 50),  # very short test runs
  static = TRUE  # with static partitioner
)               

# for additional data munging please refer to the appendix

## ---- munge-chunking-scaling, include=FALSE---------------------------------------------
scaling_chunking <- merge(scaling_chunking, chunking_bench, by = "grainsize")

single_chunk  <- transform(
  subset(scaling_chunking, chunks == 1),
  runtime_single = runtime, runtime = NULL, 
  grainsize = NULL, chunks=NULL
)

scaling_chunking <- transform(
  merge(scaling_chunking, single_chunk),
  slowdown = runtime/runtime_single,
  iter = factor(iter),
  runtime_single = NULL
)

## ---------------------------------------------------------------------------------------
ggplot(scaling_chunking) +
  aes(chunks, slowdown, colour = iter, shape = iter) +
  geom_line() + geom_point() +
  scale_x_log10(breaks = scaling_chunking$chunks) +
  scale_y_log10() +
  ggtitle("Slowdown with increasing number of chunks")

## ---- speedup-scale, message=FALSE, warning=FALSE, results='hide'-----------------------
num_cpu <- parallel::detectCores(logical = FALSE)
num_cpu_logical <- parallel::detectCores(logical = TRUE)
grainsize_default <- ceiling(N / (2 * num_cpu))
cores <- c(2^seq(0, floor(log2(num_cpu_logical))), num_cpu, num_cpu_logical)
cores <- sort(unique(cores))
grainsize <- c(2*grainsize_default, grainsize_default, grainsize_default/2)
grainsize <- round(grainsize)

scaling_cores <- benchmark_threading(
  model_poisson,
  cores = cores,
  grainsize = grainsize,
  iter = c(25), 
  static = FALSE
)

single_core  <- transform(
  subset(scaling_cores, cores == 1),
  runtime_single = runtime, runtime = NULL, cores = NULL
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
  scale_y_log10() +
  theme(legend.position = c(0.85, 0.8)) +
  ggtitle("Runtime with varying number of cores")

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
#    iter = 100, # short sampling to speedup example
#    prior = prior(normal(0,1), class = b) +
#      prior(constant(1), class = sd, group = g),
#    backend = "cmdstanr",
#    threads = threading(2)
#  )

## ---- eval=FALSE------------------------------------------------------------------------
#  # Benchmarks given model with cross-product of tuning parameters CPU
#  # cores, grainsize and iterations. Models are run with either static
#  # or non-static scheduler and inits is set by default to 0 on the
#  # unconstrained scale. Function returns a data-frame with the
#  # cross-product of the tuning parameters and as result column the
#  # respective runtime.
#  benchmark_threading <- function(model, cores = 1, grainsize = 1, iter = 100,
#                                  static = FALSE, inits = 0) {
#  
#    scaling_model <- update(
#      model, refresh = 0,
#      threads = threading(1, grainsize = grainsize[1], static = static),
#      chains = 1, iter = 2, backend = "cmdstanr"
#    )
#  
#    run_benchmark <- function(cores, size, iter) {
#      unname(
#        system.time(
#          update(
#            scaling_model, iter = iter,
#            chains = 1, seed = 1234, inits = inits, refresh = 0,
#            threads = threading(cores, grainsize = size, static = static)
#          )
#        )["elapsed"]
#      )
#    }
#  
#    cases <- expand.grid(cores = cores, grainsize = grainsize, iter = iter)
#    cases$runtime <- with(cases, mapply(run_benchmark, cores, grainsize, iter))
#    cases
#  }

## ---- eval=FALSE------------------------------------------------------------------------
#  scaling_chunking <- merge(scaling_chunking, chunking_bench, by = "grainsize")
#  
#  single_chunk  <- transform(
#    subset(scaling_chunking, chunks == 1),
#    runtime_single = runtime, runtime = NULL,
#    grainsize = NULL, chunks=NULL
#  )
#  
#  scaling_chunking <- transform(
#    merge(scaling_chunking, single_chunk),
#    slowdown = runtime/runtime_single,
#    iter = factor(iter),
#    runtime_single = NULL
#  )

