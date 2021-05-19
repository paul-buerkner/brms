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
#   normalized: suffixes of Stan lpdfs or lpmfs which only exist as normalized 
#     versions; can also be "" in which case the family is always normalized
#   specials: character vector specialties of some families

.family_gaussian <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma"), type = "real", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "se", "cens", "trunc", "mi"),
    normalized = c("_time_hom", "_time_het", "_lagsar", "_errorsar", "_fcor"),
    specials = c("residuals", "rescor")
  )
}

.family_student <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma", "nu"), type = "real", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "se", "cens", "trunc", "mi"),
    include = "fun_logm1.stan",
    normalized = c("_time_hom", "_time_het", "_lagsar", "_errorsar", "_fcor"),
    specials = c("residuals", "rescor")
  )
}

.family_skew_normal <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma", "alpha"), type = "real", 
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "se", "cens", "trunc", "mi")
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
    specials = c("categorical", "joint_link")
  )
}

.family_multinomial <- function() {
  list(
    links = "logit", 
    dpars = NULL,  # is determind based on the data
    type = "int", ybounds = c(-Inf, Inf), 
    closed = c(NA, NA),
    ad = c("weights", "subset", "trials"), 
    specials = c("multinomial", "joint_link"),
    include = "fun_multinomial_logit.stan",
    normalized = ""
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
    dpars = "phi",  # more dpars are determined based on the data
    type = "real", ybounds = c(0, 1), 
    closed = c(FALSE, FALSE),
    ad = c("weights", "subset"), 
    specials = c("dirichlet", "joint_link"),
    include = "fun_dirichlet_logit.stan",
    normalized = ""
  )
}

.family_dirichlet2 <- function() {
  list(
    links = c("log", "softplus", "squareplus", "identity", "logm1"), 
    dpars = NULL,  # is determind based on the data
    type = "real", ybounds = c(0, 1), 
    closed = c(FALSE, FALSE),
    ad = c("weights", "subset"), 
    specials = c("dirichlet"),
    include = "fun_logm1.stan",
    normalized = ""
  )
}

.family_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc", "rate")
  )
}

.family_negbinomial <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "shape"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc", "rate")
  )
}

# as negbinomial but with sigma = 1 / shape parameterization
.family_negbinomial2 <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "sigma"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc", "rate")
  )
}
 
.family_geometric <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc", "rate")
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
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "shape"), type = "int", 
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_com_poisson.stan"
  )
}

.family_gamma <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "shape"), type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_weibull <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "shape"), type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_exponential <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus", "squareplus"),
    dpars = "mu", type = "real", 
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_frechet <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "nu"), type = "real",
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = "fun_logm1.stan",
    specials = "transeta"  # see stan_eta_ilink()
  )
}

.family_inverse.gaussian <- function() {
  list(
    links = c("1/mu^2", "inverse", "identity", "log", "softplus", "squareplus"),
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
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma", "beta"), type = "real",
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi")
  )
}

.family_wiener <- function() {
  list(
    links = c("identity", "log", "softplus", "squareplus"),
    dpars = c("mu", "bs", "ndt", "bias"), type = "real",
    ybounds = c(0, Inf), closed = c(FALSE, NA),
    ad = c("weights", "subset", "dec"),
    include = "fun_wiener_diffusion.stan",
    normalized = ""
  )
}

.family_gen_extreme_value <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma", "xi"), type = "real",
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = c("fun_gen_extreme_value.stan", "fun_scale_xi.stan"),
    normalized = ""
  )
}

.family_von_mises <- function() {
  list(
    links = c("tan_half", "identity"),
    dpars = c("mu", "kappa"), type = "real",
    ybounds = c(-pi, pi), closed = c(TRUE, TRUE),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = c("fun_tan_half.stan", "fun_von_mises.stan"),
    normalized = ""
  )
}

.family_asym_laplace <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma", "quantile"), type = "real",
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc", "mi"),
    include = "fun_asym_laplace.stan",
    normalized = ""
  )
}

.family_zero_inflated_asym_laplace <- function() {
  list(
    links = c("identity", "log", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "sigma", "quantile", "zi"), type = "real",
    ybounds = c(-Inf, Inf), closed = c(NA, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = c("fun_asym_laplace.stan", "fun_zero_inflated_asym_laplace.stan")
  )
}

.family_cox <- function() {
  list(
    links = c("log", "identity", "softplus", "squareplus"),
    dpars = c("mu"), type = "real",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_cox.stan",
    specials = c("cox"),
    normalized = ""
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
    ad = c("weights", "subset", "thres", "cat"), 
    specials = c("ordinal", "ordered_thres", "thres_minus_eta", "joint_link"),
    normalized = ""
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
    ad = c("weights", "subset", "thres", "cat"), 
    specials = c("ordinal", "cs", "thres_minus_eta", "joint_link"),
    normalized = ""
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
    ad = c("weights", "subset", "thres", "cat"),
    specials = c("ordinal", "cs", "eta_minus_thres", "joint_link"),
    normalized = ""
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
    ad = c("weights", "subset", "thres", "cat"), 
    specials = c("ordinal", "cs", "eta_minus_thres", "joint_link"),
    normalized = ""
  )
}

.family_hurdle_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "hu"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_poisson.stan",
    normalized = ""
  )
}

.family_hurdle_negbinomial <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "shape", "hu"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_negbinomial.stan",
    normalized = ""
  )
}

.family_hurdle_gamma <- function() {
  list(
    links = c("log", "identity", "inverse", "softplus", "squareplus"),
    dpars = c("mu", "shape", "hu"), type = "real",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_gamma.stan",
    normalized = ""
  )
}

.family_hurdle_lognormal <- function() {
  list(
    links = c("identity", "inverse"),
    dpars = c("mu", "sigma", "hu"), type = "real",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_hurdle_lognormal.stan",
    specials = "logscale",
    normalized = ""
  )
}

.family_zero_inflated_poisson <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "zi"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_zero_inflated_poisson.stan",
    normalized = ""
  )
}

.family_zero_inflated_negbinomial <- function() {
  list(
    links = c("log", "identity", "sqrt", "softplus", "squareplus"),
    dpars = c("mu", "shape", "zi"), type = "int",
    ybounds = c(0, Inf), closed = c(TRUE, NA),
    ad = c("weights", "subset", "cens", "trunc"),
    include = "fun_zero_inflated_negbinomial.stan",
    normalized = ""
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
    include = "fun_zero_inflated_binomial.stan",
    normalized = ""
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
    include = "fun_zero_inflated_beta.stan",
    normalized = ""
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
    include = "fun_zero_one_inflated_beta.stan",
    normalized = ""
  )
}

.family_custom <- function() {
  list(
    ad = c("weights", "subset", "se", "cens", "trunc", "trials", 
           "thres", "cat", "dec", "mi", "vreal", "vint"),
    ybounds = c(-Inf, Inf), closed = c(NA, NA)
  )
}
