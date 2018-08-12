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
library(brms)
theme_set(theme_default())

## ---- include = FALSE-------------------------------------------------------------------
my_path <- "C:/Users/paulb/Dropbox/brms/Models/MCMCglmm/"

## ---------------------------------------------------------------------------------------
# uncomment and set your path to the folder where you stored the data
# my_path <- "insert your path here"
phylo <- ape::read.nexus(paste0(my_path, "phylo.nex"))
data_simple <- read.table(paste0(my_path, "data_simple.txt"), header = TRUE)
head(data_simple)

## ---------------------------------------------------------------------------------------
inv.phylo <- MCMCglmm::inverseA(phylo, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

## ---- results='hide'--------------------------------------------------------------------
model_simple <- brm(
  phen ~ cofactor + (1|phylo), data = data_simple, 
  family = gaussian(), cov_ranef = list(phylo = A),
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
plot(marginal_effects(model_simple), points = TRUE) 

## ---------------------------------------------------------------------------------------
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(model_simple, hyp, class = NULL))
plot(hyp)

## ---------------------------------------------------------------------------------------
data_repeat <- read.table(paste0(my_path, "data_repeat.txt"), header = TRUE)
data_repeat$spec_mean_cf <- 
  with(data_repeat, sapply(split(cofactor, phylo), mean)[phylo])
head(data_repeat)

## ---- results='hide'--------------------------------------------------------------------
model_repeat1 <- brm(
  phen ~ spec_mean_cf + (1|phylo) + (1|species), 
  data = data_repeat, family = gaussian(), 
  cov_ranef = list(phylo = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
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
  model_repeat1, formula = ~ . + within_spec_cf,
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
data_fisher <- read.table(paste0(my_path, "data_effect.txt"), header = TRUE)
data_fisher$obs <- 1:nrow(data_fisher)
head(data_fisher)

## ---- results='hide'--------------------------------------------------------------------
model_fisher <- brm(
  Zr | se(sqrt(1 / (N - 3))) ~ 1 + (1|phylo) + (1|obs), 
  data = data_fisher, family = gaussian(), 
  cov_ranef = list(phylo = A),
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
data_pois <- read.table(paste0(my_path, "data_pois.txt"), header = TRUE)
data_pois$obs <- 1:nrow(data_pois)
head(data_pois)

## ---- results='hide'--------------------------------------------------------------------
model_pois <- brm(
  phen_pois ~ cofactor + (1|phylo) + (1|obs), 
  data = data_pois, family = poisson("log"), 
  cov_ranef = list(phylo = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

## ---------------------------------------------------------------------------------------
summary(model_pois)
plot(marginal_effects(model_pois), points = TRUE) 

## ---- results='hide'--------------------------------------------------------------------
model_normal <- brm(
  phen_pois ~ cofactor + (1|phylo), 
  data = data_pois, family = gaussian(), 
  cov_ranef = list(phylo = A),
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

