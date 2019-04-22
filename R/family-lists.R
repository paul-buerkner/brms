# This file contains a list for every native family.
# These lists may contain the following elements:
#   links: possible link function (first is default)
#   dpars: distributional parameters of the family
#   type: either real or int (i.e. continuous or discrete)
#   ybounds: area of definition of the response values
#   closed: is the interval closed or open?
#   ad: supported addition arguments
#   include: names of user-defined Stan functions
#     to be included in the Stan code
#   const: definitions of constants in Stan
#     to be put in the transformed parameters block
#   specials: character vector specialities of some families

.family_gaussian <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus"),
    dpars = c("mu", "sigma"), type = "real", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "se", "cens", "trunc", "mi"),
    specials = c("residuals", "rescor")
  )
}

.family_student <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus"),
    dpars = c("mu", "sigma", "nu"), type = "real", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "se", "cens", "trunc", "mi"),
    include = "fun_logm1.stan",
    specials = c("residuals", "rescor")
  )
}

.family_skew_normal <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus"),
    dpars = c("mu", "sigma", "alpha"), type = "real", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "se", "cens", "trunc", "mi"),
    const = "real sqrt_2_div_pi = sqrt(2 / pi())"
  )
}

.family_binomial <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "trials", "cens", "trunc")
  )
}

.family_bernoulli <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu"), type = "int", 
    ybounds = c(0, 1), closed = c(TRUE, TRUE),
    ad = c("weights", "subset"), specials = "binary"
  )
}

.family_categorical <- function() {
  list(
    links = "logit", 
    dpars = NULL,  # is determind based on the data
    type = "int", ybounds = c(-Inf, Inf), 
    closed = c(NA, NA),
    ad = c("weights", "subset"), 
    specials = c("categorical")
  )
}

.family_multinomial <- function() {
  list(
    links = "logit", 
    dpars = NULL,  # is determind based on the data
    type = "int", ybounds = c(-Inf, Inf), 
    closed = c(NA, NA),
    ad = c("weights", "subset", "trials"), 
    specials = c("multinomial"),
    include = "fun_multinomial_logit.stan"
  )
}

.family_beta <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx",
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu", "phi"), type = "real",
    ybounds = c(0, 1), closed = c(FALSE, FALSE),
    ad = c("weights", "subset", "cens", "trunc", "mi")
  )
}

.family_dirichlet <- function() {
  list(
    links = "logit", 
    dpars = "phi",  # more dpars are determind based on the data
    type = "real", ybounds = c(0, 1), 
    closed = c(FALSE, FALSE),
    ad = c("weights", "subset"), 
    specials = c("dirichlet"),
    include = "fun_dirichlet_logit.stan"
  )
}

.family_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc")
  )
}

.family_negbinomial <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu", "shape"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc")
  )
}
 
.family_geometric <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc")
  )
}

.family_discrete_weibull <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx",
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu", "shape"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_discrete_weibull.stan"
  )
}

.family_com_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu", "shape"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_com_poisson.stan"
  )
}

.family_gamma <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus"),
    dpars = c("mu", "shape"), type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_weibull <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus"),
    dpars = c("mu", "shape"), type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_exponential <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus"),
    dpars = "mu", type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_frechet <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus"),
    dpars = c("mu", "nu"), type = "real",
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = "fun_logm1.stan",
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_inverse.gaussian <- function() {
  list(
    links = c("1/mu^2", "inverse", "identity", "log", "softplus"),
    dpars = c("mu", "shape"), type = "real",
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = "fun_inv_gaussian.stan"
  )
}

.family_lognormal <- function() {
  list(
    links = c("identity", "inverse"),
    dpars = c("mu", "sigma"), type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "logscale"
  )
}

.family_shifted_lognormal <- function() {
  list(
    links = c("identity", "inverse"),
    dpars = c("mu", "sigma", "ndt"), type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "logscale"
  )
}

.family_exgaussian <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus"),
    dpars = c("mu", "sigma", "beta"), type = "real",
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi")
  )
}

.family_wiener <- function() {
  list(
    links = "identity",
    dpars = c("mu", "bs", "ndt", "bias"), type = "real",
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "dec"),
    include = "fun_wiener_diffusion.stan"
  )
}

.family_gen_extreme_value <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus"),
    dpars = c("mu", "sigma", "xi"), type = "real",
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = c("fun_gen_extreme_value.stan", "fun_scale_xi.stan")
  )
}

.family_von_mises <- function() {
  list(
    links = "tan_half",
    dpars = c("mu", "kappa"), type = "real",
    ybounds = c(-pi, pi), closed = c(TRUE, TRUE),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = c("fun_tan_half.stan", "fun_von_mises.stan")
  )
}

.family_asym_laplace <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus"),
    dpars = c("mu", "sigma", "quantile"),
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = "fun_asym_laplace.stan"
  )
}

.family_cumulative <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit"
    ),
    dpars = c("mu", "disc"), type = "int", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cat"), 
    specials = c("ordinal", "ordered_thres", "thres_minus_eta")
  )
}

.family_sratio <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit"
    ),
    dpars = c("mu", "disc"), type = "int", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cat"), 
    specials = c("ordinal", "cs", "thres_minus_eta")
  )
}

.family_cratio <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit"
    ),
    dpars = c("mu", "disc"), type = "int", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cat"), 
    specials = c("ordinal", "cs", "eta_minus_thres")
  )
}

.family_acat <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit"
    ),
    dpars = c("mu", "disc"), type = "int", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cat"), 
    specials = c("ordinal", "cs", "eta_minus_thres")
  )
}

.family_hurdle_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu", "hu"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_poisson.stan"
  )
}

.family_hurdle_negbinomial <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu", "shape", "hu"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_negbinomial.stan"
  )
}

.family_hurdle_gamma <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus"),
    dpars = c("mu", "shape", "hu"), type = "real",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_gamma.stan"
  )
}

.family_hurdle_lognormal <- function() {
  list(
    links = c("identity", "inverse"),
    dpars = c("mu", "sigma", "hu"), type = "real",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_lognormal.stan",
    specials = "logscale"
  )
}

.family_zero_inflated_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu", "zi"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_zero_inflated_poisson.stan"
  )
}

.family_zero_inflated_negbinomial <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus"),
    dpars = c("mu", "shape", "zi"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_zero_inflated_negbinomial.stan"
  )
}

.family_zero_inflated_binomial <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx",
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu", "zi"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "trials", "cens", "trunc"),
    include = "fun_zero_inflated_binomial.stan"
  )
}

.family_zero_inflated_beta <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx",
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu", "phi", "zi"), type = "real",
    ybounds = c(0, 1), closed = c(TRUE, FALSE),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_zero_inflated_beta.stan"
  )
}

.family_zero_one_inflated_beta <- function() {
  list(
    links = c(
      "logit", "probit", "probit_approx",
      "cloglog", "cauchit", "identity"
    ),
    dpars = c("mu", "phi", "zoi", "coi"), type = "real",
    ybounds = c(0, 1), closed = c(TRUE, TRUE),
    ad = c("weights", "subset"),
    include = "fun_zero_one_inflated_beta.stan"
  )
}

.family_custom <- function() {
  list(
    ad = c("weights", "subset", "se", "cens", "trunc", "trials", "cat", "dec", "mi"),
    ybounds = c(-Inf, Inf), closed = c(NA, NA)
  )
}
