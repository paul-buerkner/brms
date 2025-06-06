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
library(brms)
ggplot2::theme_set(theme_default())

## ---------------------------------------------------------------------------------------
phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
data_simple <- read.table(
  "https://paul-buerkner.github.io/data/data_simple.txt",
  header = TRUE
)
head(data_simple)

## ---------------------------------------------------------------------------------------
A <- ape::vcv.phylo(phylo)

## ---- results='hide'--------------------------------------------------------------------
model_simple <- brm(
  phen ~ cofactor + (1 | gr(phylo, cov = A)),
  data = data_simple,
  family = gaussian(),
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  )
)

## ---------------------------------------------------------------------------------------
summary(model_simple)
plot(model_simple, N = 2, ask = FALSE)
plot(conditional_effects(model_simple), points = TRUE)

## ---------------------------------------------------------------------------------------
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(model_simple, hyp, class = NULL))
plot(hyp)

## ---------------------------------------------------------------------------------------
data_repeat <- read.table(
  "https://paul-buerkner.github.io/data/data_repeat.txt",
  header = TRUE
)
data_repeat$spec_mean_cf <-
  with(data_repeat, sapply(split(cofactor, phylo), mean)[phylo])
head(data_repeat)

## ---- results='hide'--------------------------------------------------------------------
model_repeat1 <- brm(
  phen ~ spec_mean_cf + (1 | gr(phylo, cov = A)) + (1 | species),
  data = data_repeat,
  family = gaussian(),
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

## ---------------------------------------------------------------------------------------
summary(model_repeat1)

## ---------------------------------------------------------------------------------------
hyp <- paste(
  "sd_phylo__Intercept^2 /",
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(model_repeat1, hyp, class = NULL))
plot(hyp)

## ---------------------------------------------------------------------------------------
data_repeat$within_spec_cf <- data_repeat$cofactor - data_repeat$spec_mean_cf

## ---- results='hide'--------------------------------------------------------------------
model_repeat2 <- update(
  model_repeat1,
  formula = ~ . + within_spec_cf,
  newdata = data_repeat, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

## ---------------------------------------------------------------------------------------
summary(model_repeat2)

## ---------------------------------------------------------------------------------------
hyp <- paste(
  "sd_phylo__Intercept^2 /",
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(model_repeat2, hyp, class = NULL))

## ---------------------------------------------------------------------------------------
data_fisher <- read.table(
  "https://paul-buerkner.github.io/data/data_effect.txt",
  header = TRUE
)
data_fisher$obs <- 1:nrow(data_fisher)
head(data_fisher)

## ---- results='hide'--------------------------------------------------------------------
model_fisher <- brm(
  Zr | se(sqrt(1 / (N - 3))) ~ 1 + (1 | gr(phylo, cov = A)) + (1 | obs),
  data = data_fisher, family = gaussian(),
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 10), "sd")
  ),
  control = list(adapt_delta = 0.95),
  chains = 2, cores = 2, iter = 4000, warmup = 1000
)

## ---------------------------------------------------------------------------------------
summary(model_fisher)
plot(model_fisher)

## ---------------------------------------------------------------------------------------
data_pois <- read.table(
  "https://paul-buerkner.github.io/data/data_pois.txt",
  header = TRUE
)
data_pois$obs <- 1:nrow(data_pois)
head(data_pois)

## ---- results='hide'--------------------------------------------------------------------
model_pois <- brm(
  phen_pois ~ cofactor + (1 | gr(phylo, cov = A)) + (1 | obs),
  data = data_pois, family = poisson("log"),
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

## ---------------------------------------------------------------------------------------
summary(model_pois)
plot(conditional_effects(model_pois), points = TRUE)

## ---- results='hide'--------------------------------------------------------------------
model_normal <- brm(
  phen_pois ~ cofactor + (1 | gr(phylo, cov = A)),
  data = data_pois, family = gaussian(),
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

## ---------------------------------------------------------------------------------------
summary(model_normal)

## ---------------------------------------------------------------------------------------
pp_check(model_pois)
pp_check(model_normal)

## ---------------------------------------------------------------------------------------
loo(model_pois, model_normal)
