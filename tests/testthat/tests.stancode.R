context("Tests for stancode")

# simplifies manual calling of tests
expect_match2 <- brms:::expect_match2
SW <- brms:::SW

# parsing the Stan code ensures syntactical correctness of models
# setting this option to FALSE speeds up testing
not_cran <- identical(Sys.getenv("NOT_CRAN"), "true")
options(brms.parse_stancode = not_cran, brms.backend = "rstan")

test_that("specified priors appear in the Stan code", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10),
                    g = rep(1:5, 2), h = factor(rep(1:5, each = 2)))

  prior <- c(prior(std_normal(), coef = x1),
             prior(normal(0,2), coef = x2),
             prior(normal(0,5), Intercept, lb = 0),
             prior(cauchy(0,1), sd, group = g, lb = "", ub = 5),
             prior(cauchy(0,2), sd, group = g, coef = x1),
             prior(gamma(1, 1), sd, group = h, ub = 10))
  scode <- stancode(y ~ x1*x2 + (x1*x2|g) + (1 | h), dat,
                         prior = prior, sample_prior = "yes")
  expect_match2(scode, "vector<upper=5>[M_1] sd_1;")
  expect_match2(scode, "vector<lower=0,upper=10>[M_2] sd_2;")
  expect_match2(scode, "target += lprior;")
  expect_match2(scode, "lprior += std_normal_lpdf(b[1])")
  expect_match2(scode, "lprior += normal_lpdf(b[2] | 0, 2)")
  expect_match2(scode, "lprior += normal_lpdf(Intercept | 0, 5)")
  expect_match2(scode, "lprior += cauchy_lpdf(sd_1[1] | 0, 1)")
  expect_match2(scode, "- 1 * cauchy_lcdf(5 | 0, 1)")
  expect_match2(scode, "lprior += cauchy_lpdf(sd_1[2] | 0, 2)")
  expect_match2(scode, "lprior += student_t_lpdf(sigma | 3, 0, 3.7)")
  expect_match2(scode, "- 1 * student_t_lccdf(0 | 3, 0, 3.7)")
  expect_match2(scode, "lprior += gamma_lpdf(sd_2 | 1, 1)")
  expect_match2(scode, "prior_b__1 = normal_rng(0,1);")
  expect_match2(scode, "prior_sd_1__1 = cauchy_rng(0,1)")
  expect_match2(scode, "while (prior_sd_1__1 > 5)")
  expect_match2(scode, "prior_sd_2 = gamma_rng(1,1)")
  expect_match2(scode, "while (prior_sd_2 < 0 || prior_sd_2 > 10)")

  prior <- c(prior(lkj(0.5), class = cor, group = g),
             prior(normal(0, 1), class = b),
             prior(normal(0, 5), class = Intercept),
             prior(cauchy(0, 5), class = sd))
  scode <- stancode(y ~ x1 + cs(x2) + (0 + x1 + x2 | g),
                         data = dat, family = acat(),
                         prior = prior, sample_prior = TRUE)
  expect_match2(scode, "lprior += normal_lpdf(b | 0, 1)")
  expect_match2(scode, "lprior += normal_lpdf(Intercept | 0, 5)")
  expect_match2(scode, "lprior += cauchy_lpdf(sd_1 | 0, 5)")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(L_1 | 0.5)")
  expect_match2(scode, "lprior += normal_lpdf(to_vector(bcs) | 0, 1)")
  expect_match2(scode, "prior_bcs = normal_rng(0,1)")

  prior <- c(prior(normal(0,5), nlpar = a),
             prior(normal(0,10), nlpar = b),
             prior(cauchy(0,1), class = sd, nlpar = a),
             prior(lkj(2), class = cor, group = g))
  scode <- stancode(
    bf(y ~ a * exp(-b * x1), a + b ~ (1|ID|g), nl = TRUE),
    data = dat, prior = prior, sample_prior = TRUE
  )
  expect_match2(scode, "lprior += normal_lpdf(b_a | 0, 5)")
  expect_match2(scode, "lprior += normal_lpdf(b_b | 0, 10)")
  expect_match2(scode, "lprior += cauchy_lpdf(sd_1[1] | 0, 1)")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(L_1 | 2)")
  expect_match2(scode, "prior_b_a = normal_rng(0,5)")
  expect_match2(scode, "prior_sd_1__2 = student_t_rng(3,0,3.7)")
  expect_match2(scode, "prior_cor_1 = lkj_corr_rng(M_1,2)[1, 2]")

  prior <- c(prior(lkj(2), rescor),
             prior(cauchy(0, 5), sigma, resp = y),
             prior(cauchy(0, 1), sigma, resp = x1))
  form <- bf(mvbind(y, x1) ~ x2) + set_rescor(TRUE)
  scode <- stancode(form, dat, prior = prior, sample_prior = TRUE)
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(Lrescor | 2)")
  expect_match2(scode, "prior_sigma_y = cauchy_rng(0,5)")
  expect_match2(scode, "prior_rescor = lkj_corr_rng(nresp,2)[1, 2]")

  prior <- c(prior(uniform(-1, 1), ar),
             prior(normal(0, 0.5), ma),
             prior(normal(0, 5)))
  scode <- stancode(y ~ mo(g) + arma(cov = TRUE), dat,
                         prior = prior, sample_prior = TRUE)
  expect_match2(scode, "vector<lower=-1,upper=1>[Kar] ar;")
  expect_match2(scode, "vector<lower=-1,upper=1>[Kma] ma;")
  expect_match2(scode, "lprior += uniform_lpdf(ar | -1, 1)")
  expect_match2(scode, "lprior += normal_lpdf(ma | 0, 0.5)")
  expect_match2(scode,
    "- 1 * log_diff_exp(normal_lcdf(1 | 0, 0.5), normal_lcdf(-1 | 0, 0.5))"
  )
  expect_match2(scode, "lprior += normal_lpdf(bsp | 0, 5)")
  expect_match2(scode, "lprior += dirichlet_lpdf(simo_1 | con_simo_1)")
  expect_match2(scode, "prior_simo_1 = dirichlet_rng(con_simo_1)")
  expect_match2(scode, "prior_ar = uniform_rng(-1,1)")
  expect_match2(scode, "while (prior_ar < -1 || prior_ar > 1)")

  # test for problem described in #213
  prior <- c(prior(normal(0, 1), coef = x1),
             prior(normal(0, 2), coef = x1, dpar = sigma))
  scode <- stancode(bf(y ~ x1, sigma ~ x1), dat, prior = prior)
  expect_match2(scode, "lprior += normal_lpdf(b[1] | 0, 1);")
  expect_match2(scode, "lprior += normal_lpdf(b_sigma[1] | 0, 2);")

  prior <- c(set_prior("target += normal_lpdf(b[1] | 0, 1)", check = FALSE),
             set_prior("", class = "sigma"))
  scode <- stancode(y ~ x1, dat, prior = prior, sample_prior = TRUE)
  expect_match2(scode, "target += normal_lpdf(b[1] | 0, 1)")
  expect_true(!grepl("sigma \\|", scode))

  # tests use of default priors stored in families #1614
  scode <- stancode(y ~ x1, dat, family = negbinomial())
  expect_match2(scode, "lprior += inv_gamma_lpdf(shape | 0.4, 0.3);")
  scode <- stancode(y ~ x2, dat, family = von_mises())
  expect_match2(scode, "lprior += student_t_lpdf(Intercept | 1, 0, 1);")

  prior <- prior(gamma(0, 1), coef = x1)
  expect_warning(stancode(y ~ x1, dat, prior = prior),
                 "no natural lower bound")
  prior <- prior(uniform(0,5), class = sd)
  expect_warning(stancode(y ~ x1 + (1|g), dat, prior = prior),
                  "no natural upper bound")

  prior <- prior(uniform(-1, 1), class = cor)
  expect_error(
    stancode(y ~ x1 + (x1|g), dat, prior = prior),
    "prior for correlation matrices is the 'lkj' prior"
  )
})

test_that("special shrinkage priors appear in the Stan code", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10),
                    g = rep(1:2, each = 5), x3 = sample(1:5, 10, TRUE))

  # horseshoe prior
  hs <- horseshoe(7, scale_global = 2, df_global = 3,
                  df_slab = 6, scale_slab = 3)
  scode <- stancode(y ~ x1*x2, data = dat,
                         prior = set_prior(hs),
                         sample_prior = TRUE)
  expect_match2(scode, "vector<lower=0>[Kscales] hs_local;")
  expect_match2(scode, "real<lower=0> hs_global;")
  expect_match2(scode,
    "target += student_t_lpdf(hs_local | hs_df, 0, 1)"
  )
  expect_match2(scode,
    "lprior += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global * sigma)"
  )
  expect_match2(scode,
    "lprior += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab)"
  )
  expect_match2(scode,
    "scales = scales_horseshoe(hs_local, hs_global, hs_scale_slab^2 * hs_slab);"
  )

  scode <- stancode(y ~ x1*x2, data = dat, poisson(),
                         prior = prior(horseshoe(scale_global = 3)))
  expect_match2(scode,
    "scales = scales_horseshoe(hs_local, hs_global, hs_scale_slab^2 * hs_slab);"
  )

  scode <- stancode(x1 ~ mo(y), dat, prior = prior(horseshoe()))
  expect_match2(scode, "target += std_normal_lpdf(zbsp);")
  expect_match2(scode,
    "target += student_t_lpdf(hs_local | hs_df, 0, 1)"
  )
  expect_match2(scode,
    "scales = scales_horseshoe(hs_local, hs_global, hs_scale_slab^2 * hs_slab);"
  )

  # R2D2 prior
  scode <- stancode(y ~ x1*x2, data = dat,
                         prior = prior(R2D2(0.5, 10)),
                         sample_prior = TRUE)
  expect_match2(scode, "scales = scales_R2D2(R2D2_phi, R2D2_tau2);")
  expect_match2(scode, "target += dirichlet_lpdf(R2D2_phi | R2D2_cons_D2);")
  expect_match2(scode, "lprior += beta_lpdf(R2D2_R2 | R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2);")
  expect_match2(scode, "R2D2_tau2 = sigma^2 * R2D2_R2 / (1 - R2D2_R2);")

  # shrinkage priors applied in a non-linear model
  hs_a1 <- horseshoe(7, scale_global = 2, df_global = 3)
  R2D2_a2 <- R2D2(0.5, 10)
  scode <- SW(stancode(
    bf(y ~ a1 + a2, a1 ~ x1, a2 ~ 0 + x2, nl = TRUE),
    data = dat, sample_prior = TRUE,
    prior = c(set_prior(hs_a1, nlpar = "a1"),
              set_prior(R2D2_a2, nlpar = "a2"))
  ))
  expect_match2(scode, "vector<lower=0>[Kscales_a1] hs_local_a1;")
  expect_match2(scode, "real<lower=0> hs_global_a1;")
  expect_match2(scode,
    "target += student_t_lpdf(hs_local_a1 | hs_df_a1, 0, 1)"
  )
  expect_match2(scode,
    "lprior += student_t_lpdf(hs_global_a1 | hs_df_global_a1, 0, hs_scale_global_a1 * sigma)"
  )
  expect_match2(scode,
    "lprior += inv_gamma_lpdf(hs_slab_a1 | 0.5 * hs_df_slab_a1, 0.5 * hs_df_slab_a1)"
  )
  expect_match2(scode,
    "scales_a1 = scales_horseshoe(hs_local_a1, hs_global_a1, hs_scale_slab_a1^2 * hs_slab_a1);"
  )
  expect_match2(scode, "scales_a2 = scales_R2D2(R2D2_phi_a2, R2D2_tau2_a2);")

  # shrinkage priors can be applied globally
  bform <- bf(y ~ x1*mo(x3) + (1|g) + (1|x1) + gp(x3) + s(x2) +
                arma(p = 2, q = 2, gr = g))
  bprior <- prior(R2D2(main = TRUE), class = b) +
    prior(R2D2(), class = sd) +
    prior(R2D2(), class = sds) +
    prior(R2D2(), class = sdgp) +
    prior(R2D2(), class = ar) +
    prior(R2D2(), class = ma)
  scode <- stancode(bform, data = dat, prior = bprior)
  expect_match2(scode, "sdb = scales[(1):(Kc)];")
  expect_match2(scode, "sdbsp = scales[(1+Kc):(Kc+Ksp)];")
  expect_match2(scode, "sdbs = scales[(1+Kc+Ksp):(Kc+Ksp+Ks)];")
  expect_match2(scode, "sds_1 = scales[(1+Kc+Ksp+Ks):(Kc+Ksp+Ks+nb_1)];")
  expect_match2(scode, "sdgp_1 = scales[(1+Kc+Ksp+Ks+nb_1):(Kc+Ksp+Ks+nb_1+Kgp_1)];")
  expect_match2(scode, "sdar = scales[(1+Kc+Ksp+Ks+nb_1+Kgp_1):(Kc+Ksp+Ks+nb_1+Kgp_1+Kar)];")
  expect_match2(scode, "sdma = scales[(1+Kc+Ksp+Ks+nb_1+Kgp_1+Kar):(Kc+Ksp+Ks+nb_1+Kgp_1+Kar+Kma)];")
  expect_match2(scode, "sd_1 = scales[(1+Kc+Ksp+Ks+nb_1+Kgp_1+Kar+Kma):(Kc+Ksp+Ks+nb_1+Kgp_1+Kar+Kma+M_1)];")
  expect_match2(scode, "sd_2 = scales[(1+Kc+Ksp+Ks+nb_1+Kgp_1+Kar+Kma+M_1):(Kc+Ksp+Ks+nb_1+Kgp_1+Kar+Kma+M_1+M_2)];")
  expect_match2(scode, "bsp = zbsp .* sdbsp;  // scale coefficients")
  expect_match2(scode, "ar = zar .* sdar;  // scale coefficients")

  # check error messages
  expect_error(stancode(y ~ x1*x2, data = dat,
                             prior = prior(horseshoe(-1))),
               "Degrees of freedom of the local priors")
  expect_error(stancode(y ~ x1*x2, data = dat,
                             prior = prior(horseshoe(1, -1))),
               "Scale of the global prior")
  expect_error(stancode(y ~ cs(x1), dat, acat(), prior = prior(R2D2())),
               "Special priors are not yet allowed")
  bprior <- prior(horseshoe()) + prior(normal(0, 1), coef = "y")
  expect_error(stancode(x1 ~ y, dat, prior = bprior),
               "Defining separate priors for single coefficients")
  expect_error(stancode(x1 ~ y, dat, prior = prior(horseshoe(), lb = 0)),
               "Setting boundaries on coefficients is not allowed")
  expect_error(
    stancode(y ~ x1*x2, data = dat, prior = prior(lasso(2, scale = 10))),
    "The lasso prior is no longer supported"
  )
})

test_that("priors can be fixed to constants", {
  dat <- data.frame(y = 1:12, x1 = rnorm(12), x2 = rnorm(12),
                    g = rep(1:6, each = 2), h = factor(rep(1:2, each = 6)))

  prior <- prior(normal(0, 1), b) +
    prior(constant(3), b, coef = x1) +
    prior(constant(-1), b, coef = x2) +
    prior(constant(10), Intercept) +
    prior(normal(0, 5), sd) +
    prior(constant(1), sd, group = g, coef = x2) +
    prior(constant(2), sd, group = g, coef = x1) +
    prior(constant(0.3), sigma)
  scode <- stancode(y ~ x1*x2 + (x1*x2 | g), dat, prior = prior)
  expect_match2(scode, "b[1] = 3;")
  expect_match2(scode, "b[2] = -1;")
  expect_match2(scode, "b[3] = par_b_3;")
  expect_match2(scode, "lprior += normal_lpdf(b[3] | 0, 1);")
  expect_match2(scode, "Intercept = 1")
  expect_match2(scode, "sd_1[3] = 1;")
  expect_match2(scode, "sd_1[2] = 2;")
  expect_match2(scode, "sd_1[4] = par_sd_1_4;")
  expect_match2(scode, "lprior += normal_lpdf(sd_1[4] | 0, 5)")
  expect_match2(scode, "sigma = 0.3;")

  prior <- prior(constant(3))
  scode <- stancode(y ~ x2 + x1 + cs(g), dat, family = sratio(),
                         prior = prior)
  expect_match2(scode, "b = rep_vector(3, rows(b));")
  expect_match2(scode, "bcs = rep_matrix(3, rows(bcs), cols(bcs));")

  prior <- prior(normal(0, 3)) +
    prior(constant(3), coef = x1) +
    prior(constant(-1), coef = g)
  scode <- stancode(y ~ x1 + cs(x2) + cs(g), dat, family = sratio(),
                         prior = prior)
  expect_match2(scode, "b[1] = 3;")
  expect_match2(scode, "bcs[1] = par_bcs_1;")
  expect_match2(scode, "lprior += normal_lpdf(bcs[1] | 0, 3);")
  expect_match2(scode, "bcs[2] = rep_row_vector(-1, cols(bcs[2]));")

  prior <- prior(constant(3), class = "sd", group = "g") +
    prior(constant("[[1, 0], [0, 1]]"), class = "cor")
  scode <- stancode(y ~ x1 + (x1 | gr(g, by = h)), dat, prior = prior)
  expect_match2(scode, "sd_1 = rep_matrix(3, rows(sd_1), cols(sd_1));")
  expect_match2(scode, "L_1[2] = [[1, 0], [0, 1]];")

  prior <- prior(constant(0.5), class = lscale, coef = gpx1h1) +
    prior(normal(0, 10), class = lscale, coef = gpx1h2)
  scode <- stancode(y ~ gp(x1, by = h), dat, prior = prior)
  expect_match2(scode, "lscale_1[1][1] = 0.5;")
  expect_match2(scode, "lscale_1[2][1] = par_lscale_1_2_1;")
  expect_match2(scode, "lprior += normal_lpdf(lscale_1[2][1] | 0, 10)")

  # test that improper base priors are correctly recognized (#919)
  prior <- prior(constant(-1), b, coef = x2)
  scode <- stancode(y ~ x1*x2, dat, prior = prior)
  expect_match2(scode, "real par_b_1;")
  expect_match2(scode, "b[3] = par_b_3;")

  # test error messages
  prior <- prior(normal(0, 1), Intercept) +
    prior(constant(3), Intercept, coef = 2)
  expect_error(
    stancode(y ~ x1, data = dat, family = cumulative(), prior = prior),
    "Can either estimate or fix all values"
  )
})

test_that("link functions appear in the Stan code", {
  dat <- data.frame(y = 1:10, x = rnorm(10))
  expect_match2(stancode(y ~ s(x), dat, family = poisson()),
               "target += poisson_log_lpmf(Y | mu);")
  expect_match2(stancode(mvbind(y, y + 1) ~ x, dat,
                              family = skew_normal("log")),
               "mu_y = exp(mu_y);")
  expect_match2(stancode(y ~ x, dat, family = von_mises(tan_half)),
               "mu = inv_tan_half(mu);")
  expect_match2(stancode(y ~ x, dat, family = weibull()),
                "mu = exp(mu);")
  expect_match2(stancode(y ~ x, dat, family = poisson("sqrt")),
               "mu = square(mu);")
  expect_match2(stancode(y ~ s(x), dat, family = bernoulli()),
                "target += bernoulli_logit_lpmf(Y | mu);")

  scode <- stancode(y ~ x, dat, family = beta_binomial('logit'))
  expect_match2(scode, "mu = inv_logit(mu);")
  scode <- stancode(y ~ x, dat, family = beta_binomial('cloglog'))
  expect_match2(scode, "mu = inv_cloglog(mu);")
  scode <- stancode(y ~ x, dat, family = beta_binomial('cauchit'))
  expect_match2(scode, "mu = inv_cauchit(mu);")

  scode <- stancode(y ~ x, dat, family = cumulative('cauchit'))
  expect_match2(scode, "p = inv_cauchit(disc * (thres[1] - mu));")
})

test_that("Stan GLM primitives are applied correctly", {
  dat <- data.frame(x = rnorm(10), y = 1:10)

  scode <- stancode(y ~ x, dat, family = gaussian)
  expect_match2(scode, "normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma)")

  scode <- stancode(y ~ x, dat, family = bernoulli)
  expect_match2(scode, "bernoulli_logit_glm_lpmf(Y | Xc, Intercept, b)")

  scode <- stancode(y ~ x, dat, family = poisson)
  expect_match2(scode, "poisson_log_glm_lpmf(Y | Xc, Intercept, b)")

  scode <- stancode(y ~ x, dat, family = negbinomial)
  expect_match2(scode,
    "neg_binomial_2_log_glm_lpmf(Y | Xc, Intercept, b, shape)"
  )

  scode <- stancode(y ~ x, dat, family = brmsfamily("negbinomial2"))
  expect_match2(scode,
    "neg_binomial_2_log_glm_lpmf(Y | Xc, Intercept, b, inv(sigma))"
  )

  scode <- stancode(y ~ 0 + x, dat, family = gaussian)
  expect_match2(scode, "normal_id_glm_lpdf(Y | X, 0, b, sigma)")

  bform <- bf(y ~ x) + bf(x ~ 1, family = negbinomial()) + set_rescor(FALSE)
  scode <- stancode(bform, dat, family = gaussian)
  expect_match2(scode,
    "normal_id_glm_lpdf(Y_y | Xc_y, Intercept_y, b_y, sigma_y)"
  )

  scode <- stancode(bf(y ~ x, decomp = "QR"), dat, family = gaussian)
  expect_match2(scode, "normal_id_glm_lpdf(Y | XQ, Intercept, bQ, sigma);")
})

test_that("customized covariances appear in the Stan code", {
  M <- diag(1, nrow = length(unique(inhaler$subject)))
  rownames(M) <- unique(inhaler$subject)
  dat2 <- list(M = M)

  scode <- stancode(rating ~ treat + (1 | gr(subject, cov = M)),
                         data = inhaler, data2 = dat2)
  expect_match2(scode, "r_1_1 = (sd_1[1] * (Lcov_1 * z_1[1]))")

  scode <- stancode(rating ~ treat + (1 + treat | gr(subject, cov = M)),
                         data = inhaler, data2 = dat2)
  expect_match2(scode, "r_1 = scale_r_cor_cov(z_1, sd_1, L_1, Lcov_1);")
  expect_match2(scode, "cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];")

  scode <- stancode(rating ~ (1 + treat | gr(subject, cor = FALSE, cov = M)),
                         data = inhaler, data2 = dat2)
  expect_match2(scode, "r_1_1 = (sd_1[1] * (Lcov_1 * z_1[1]));")
  expect_match2(scode, "r_1_2 = (sd_1[2] * (Lcov_1 * z_1[2]));")

  inhaler$by <- inhaler$subject %% 2
  scode <- stancode(rating ~ (1 + treat | gr(subject, by = by, cov = M)),
                         data = inhaler, data2 = dat2)
  expect_match2(scode, "r_1 = scale_r_cor_by_cov(z_1, sd_1, L_1, Jby_1, Lcov_1);")

  expect_warning(
    scode <- stancode(rating ~ treat + period + carry + (1|subject),
                           data = inhaler, cov_ranef = list(subject = 1)),
    "Argument 'cov_ranef' is deprecated"
  )
  expect_match2(scode, "r_1_1 = (sd_1[1] * (Lcov_1 * z_1[1]))")
})

test_that("truncation appears in the Stan code", {
  scode <- stancode(time | trunc(0) ~ age + sex + disease,
                         data = kidney, family = "gamma")
  expect_match2(scode, "target += gamma_lpdf(Y[n] | shape, shape ./ mu[n]) -")
  expect_match2(scode, "gamma_lccdf(lb[n] | shape, shape ./ mu[n]);")

  scode <- stancode(time | trunc(ub = 100) ~ age + sex + disease,
                         data = kidney, family = student("log"))

  expect_match2(scode, "target += student_t_lpdf(Y[n] | nu, mu[n], sigma) -")
  expect_match2(scode, "student_t_lcdf(ub[n] | nu, mu[n], sigma);")

  scode <- stancode(count | trunc(0, 150) ~ Trt,
                         data = epilepsy, family = "poisson")
  expect_match2(scode, "target += poisson_lpmf(Y[n] | mu[n]) -")
  expect_match2(scode,
    "log_diff_exp(poisson_lcdf(ub[n] | mu[n]), poisson_lcdf(lb[n] - 1 | mu[n]));"
  )
})

test_that("stancode handles models without fixed effects", {
  expect_match2(stancode(count ~ 0 + (1|patient) + (1+Trt|visit),
                             data = epilepsy, family = "poisson"),
               "mu = rep_vector(0.0, N);")
})

test_that("stancode correctly restricts FE parameters", {
  data <- data.frame(y = rep(0:1, each = 5), x = rnorm(10))

  scode <- stancode(y ~ x, data, prior = set_prior("", lb = 2))
  expect_match2(scode, "vector<lower=2>[Kc] b")

  scode <- stancode(
    y ~ x, data, prior = set_prior("normal (0, 2)", ub = "4")
  )
  expect_match2(scode, "vector<upper=4>[Kc] b")
  expect_match2(scode, "- 1 * normal_lcdf(4 | 0, 2)")

  prior <- set_prior("normal(0,5)", lb = "-3", ub = 5)
  scode <- stancode(y ~ 0 + x, data, prior = prior)
  expect_match2(scode, "vector<lower=-3,upper=5>[K] b")
})

test_that("self-defined functions appear in the Stan code", {
  # cauchit link
  scode <- stancode(rating ~ treat, data = inhaler,
                         family = bernoulli("cauchit"))
  expect_match2(scode, "real inv_cauchit(real y)")

  # softplus link
  scode <- stancode(rating ~ treat, data = inhaler,
                         family = brmsfamily("poisson", "softplus"))
  expect_match2(scode, "vector log_expm1(vector x)")

  # squareplus link
  scode <- stancode(rating ~ treat, data = inhaler,
                         family = brmsfamily("poisson", "squareplus"))
  expect_match2(scode, "real squareplus(real x)")

  # tan_half link
  expect_match2(stancode(rating ~ treat, data = inhaler,
                              family = von_mises("tan_half")),
               "vector inv_tan_half(vector y)")

  # logm1 link
  expect_match2(stancode(rating ~ treat, data = inhaler,
                              family = frechet()),
                "real expp1(real y)")

  # inverse gaussian models
  scode <- stancode(time | cens(censored) ~ age, data = kidney,
                                 family = inverse.gaussian)
  expect_match2(scode, "real inv_gaussian_lpdf(real y")
  expect_match2(scode, "real inv_gaussian_lcdf(real y")
  expect_match2(scode, "real inv_gaussian_lccdf(real y")
  expect_match2(scode, "real inv_gaussian_lpdf(vector y")

  # zero-inflated and hurdle models
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                              family = "zero_inflated_poisson"),
               "real zero_inflated_poisson_lpmf(int y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                             family = "zero_inflated_negbinomial"),
               "real zero_inflated_neg_binomial_lpmf(int y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                              family = "zero_inflated_binomial"),
               "real zero_inflated_binomial_lpmf(int y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                              family = "zero_inflated_beta_binomial"),
                "real zero_inflated_beta_binomial_lpmf(int y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                             family = "zero_inflated_beta"),
               "real zero_inflated_beta_lpdf(real y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                              family = "zero_one_inflated_beta"),
                "real zero_one_inflated_beta_lpdf(real y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                             family = hurdle_poisson()),
               "real hurdle_poisson_lpmf(int y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                             family = hurdle_negbinomial),
               "real hurdle_neg_binomial_lpmf(int y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                             family = hurdle_gamma("log")),
               "real hurdle_gamma_lpdf(real y")
  expect_match2(stancode(count ~ Trt, data = epilepsy,
                             family = hurdle_lognormal("identity")),
               "real hurdle_lognormal_lpdf(real y")

  # linear models with special covariance structures
  expect_match2(
    stancode(rating ~ treat + ar(cov = TRUE), data = inhaler),
    "real normal_time_hom_lpdf(vector y"
  )
  expect_match2(
    stancode(time ~ age + ar(cov = TRUE), data = kidney,
                  family = "student"),
    "real student_t_time_hom_lpdf(vector y"
  )

  # ARMA covariance matrices
  expect_match2(
    stancode(rating ~ treat + ar(cov = TRUE), data = inhaler),
    "matrix cholesky_cor_ar1(real ar"
  )
  expect_match2(
    stancode(time ~ age + ma(cov = TRUE), data = kidney),
    "matrix cholesky_cor_ma1(real ma"
  )
  expect_match2(
    stancode(time ~ age + arma(cov = TRUE), data = kidney),
    "matrix cholesky_cor_arma1(real ar, real ma"
  )
})

test_that("invalid combinations of modeling options are detected", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10),
                     wi = 1:10, ci = sample(-1:1, 10, TRUE))
  expect_error(
    stancode(y1 | cens(ci) ~ y2 + ar(cov = TRUE), data = data),
    "Invalid addition arguments for this model"
  )
  form <- bf(mvbind(y1, y2) ~ 1 + ar(cov = TRUE)) + set_rescor(TRUE)
  expect_error(
    stancode(form, data = data),
    "Explicit covariance terms cannot be modeled when 'rescor'"
  )
  expect_error(
    stancode(y1 | resp_se(wi) ~ y2 + ma(), data = data),
    "Please set cov = TRUE in ARMA structures"
  )
})

test_that("Stan code for multivariate models is correct", {
  dat <- data.frame(
    y1 = rnorm(10), y2 = rnorm(10),
    x = 1:10, g = rep(1:2, each = 5),
    censi = sample(0:1, 10, TRUE)
  )
  # models with residual correlations
  form <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
  prior <- prior(horseshoe(2), resp = "y1") +
    prior(horseshoe(2), resp = "y2")
  scode <- stancode(form, dat, prior = prior)
  expect_match2(scode, "target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);")
  expect_match2(scode, "LSigma = diag_pre_multiply(sigma, Lrescor);")
  expect_match2(scode, "target += student_t_lpdf(hs_local_y1 | hs_df_y1, 0, 1)")
  expect_match2(scode, "target += student_t_lpdf(hs_local_y2 | hs_df_y2, 0, 1)")
  expect_match2(scode, "rescor[choose(k - 1, 2) + j] = Rescor[j, k];")

  form <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
  prior <- prior(R2D2(0.2, 10), resp = "y1") +
    prior(R2D2(0.5, 10), resp = "y2")
  scode <- SW(stancode(form, dat, student(), prior = prior))
  expect_match2(scode, "target += multi_student_t_lpdf(Y | nu, Mu, Sigma);")
  expect_match2(scode, "matrix[nresp, nresp] Sigma = multiply_lower")
  expect_match2(scode, "lprior += gamma_lpdf(nu | 2, 0.1)")
  expect_match2(scode, "target += dirichlet_lpdf(R2D2_phi_y2 | R2D2_cons_D2_y2);")

  form <- bf(mvbind(y1, y2) |  weights(x) ~ 1) + set_rescor(TRUE)
  scode <- stancode(form, dat)
  expect_match2(scode,
    "target += weights[n] * (multi_normal_cholesky_lpdf(Y[n] | Mu[n], LSigma));"
  )

  # models without residual correlations
  expect_warning(
    bform <- bf(y1 | cens(censi) ~ x + y2 + (x|2|g)) +
      gaussian() + cor_ar() +
      (bf(x ~ 1) + mixture(poisson, nmix = 2)) +
      (bf(y2 ~ s(y2) + (1|2|g)) + skew_normal()),
    "Using 'cor_brms' objects for 'autocor' is deprecated"
  )
  bprior <- prior(normal(0, 5), resp = y1) +
    prior(normal(0, 10), resp = y2)
  scode <- stancode(bform, dat, prior = bprior)
  expect_match2(scode, "r_1_y2_3 = r_1[, 3]")
  expect_match2(scode, "err_y1[n] = Y_y1[n] - mu_y1[n]")
  expect_match2(scode, "target += normal_lccdf(Y_y1[Jrcens_y1[1:Nrcens_y1]] | mu_y1[Jrcens_y1[1:Nrcens_y1]], sigma_y1)")
  expect_match2(scode, "target += skew_normal_lpdf(Y_y2 | mu_y2, omega_y2, alpha_y2)")
  expect_match2(scode, "ps[1] = log(theta1_x) + poisson_log_lpmf(Y_x[n] | mu1_x[n])")
  expect_match2(scode, "lprior += normal_lpdf(b_y1 | 0, 5)")
  expect_match2(scode, "lprior += normal_lpdf(bs_y2 | 0, 10)")

  # multivariate binomial models
  bform <- bf(x ~ 1) + bf(g ~ 1) + binomial()
  scode <- stancode(bform, dat)
  expect_match2(scode, "binomial_logit_lpmf(Y_x | trials_x, mu_x)")
  expect_match2(scode, "binomial_logit_lpmf(Y_g | trials_g, mu_g)")

  # multivariate weibull models
  bform <- bform + weibull()
  scode <- stancode(bform, dat)
  expect_match2(scode, "weibull_lpdf(Y_g | shape_g, mu_g ./ tgamma(1 + 1 ./ shape_g));")
})

test_that("Stan code for categorical models is correct", {
  dat <- data.frame(y = rep(c(1, 2, 3, "a_b"), 2), x = 1:8, .g = 1:8)
  prior <- prior(normal(0, 5), "b", dpar = muab) +
    prior(normal(0, 10), "b", dpar = mu2) +
    prior(cauchy(0, 1), "Intercept", dpar = mu2) +
    prior(normal(0, 2), "Intercept", dpar = mu3)

  scode <- stancode(y ~ x + (1 | gr(.g, id = "ID")), data = dat,
                         family = categorical(), prior = prior)
  expect_match2(scode, "target += categorical_logit_lpmf(Y[n] | mu[n]);")
  expect_match2(scode, "mu[n] = transpose([0, mu2[n], mu3[n], muab[n]]);")
  expect_match2(scode, "mu2 += Intercept_mu2 + Xc_mu2 * b_mu2;")
  expect_match2(scode, "muab[n] += r_1_muab_3[J_1[n]] * Z_1_muab_3[n];")
  expect_match2(scode, "lprior += normal_lpdf(b_mu2 | 0, 10);")
  expect_match2(scode, "lprior += normal_lpdf(b_muab | 0, 5);")
  expect_match2(scode, "lprior += cauchy_lpdf(Intercept_mu2 | 0, 1);")
  expect_match2(scode, "lprior += normal_lpdf(Intercept_mu3 | 0, 2);")
  expect_match2(scode, "r_1 = scale_r_cor(z_1, sd_1, L_1);")

  scode <- stancode(y ~ x + (1 |ID| .g), data = dat,
                         family = categorical(refcat = NA))
  expect_match2(scode, "mu[n] = transpose([mu1[n], mu2[n], mu3[n], muab[n]]);")

  # test use of glm primitive
  scode <- stancode(y ~ x, data = dat, family = categorical())
  expect_match2(scode, "b[, 1] = rep_vector(0, Kc_mu2);")
  expect_match2(scode, "b[, 3] = b_mu3;")
  expect_match2(scode, "Intercept = transpose([0, Intercept_mu2, Intercept_mu3, Intercept_muab]);")
  expect_match2(scode, "target += categorical_logit_glm_lpmf(Y | Xc_mu2, Intercept, b);")

  scode <- stancode(bf(y ~ x, center = FALSE), data = dat, family = categorical())
  expect_match2(scode, "target += categorical_logit_glm_lpmf(Y | X_mu2, rep_vector(0, ncat), b);")
})

test_that("Stan code for multinomial models is correct", {
  N <- 15
  dat <- data.frame(
    y1 = rbinom(N, 10, 0.3), y2 = rbinom(N, 10, 0.5),
    y3 = rbinom(N, 10, 0.7), x = rnorm(N)
  )
  dat$size <- with(dat, y1 + y2 + y3)
  dat$y <- with(dat, cbind(y1, y2, y3))
  prior <- prior(normal(0, 10), "b", dpar = muy2) +
    prior(cauchy(0, 1), "Intercept", dpar = muy2) +
    prior(normal(0, 2), "Intercept", dpar = muy3)
  scode <- stancode(bf(y | trials(size)  ~ 1, muy2 ~ x), data = dat,
                         family = multinomial(), prior = prior)
  expect_match2(scode, "array[N, ncat] int Y;")
  expect_match2(scode, "target += multinomial_logit2_lpmf(Y[n] | mu[n]);")
  expect_match2(scode, "muy2 += Intercept_muy2 + Xc_muy2 * b_muy2;")
  expect_match2(scode, "lprior += normal_lpdf(b_muy2 | 0, 10);")
  expect_match2(scode, "lprior += cauchy_lpdf(Intercept_muy2 | 0, 1);")
  expect_match2(scode, "lprior += normal_lpdf(Intercept_muy3 | 0, 2);")
})

test_that("Stan code for dirichlet_multinomial models is correct", {
  N <- 15
  dat <- data.frame(
    y1 = rbinom(N, 10, 0.3), y2 = rbinom(N, 10, 0.5),
    y3 = rbinom(N, 10, 0.7), x = rnorm(N)
  )
  dat$size <- with(dat, y1 + y2 + y3)
  dat$y <- with(dat, cbind(y1, y2, y3))
  prior <- prior(normal(0, 10), "b", dpar = muy2) +
    prior(cauchy(0, 1), "Intercept", dpar = muy2) +
    prior(normal(0, 2), "Intercept", dpar = muy3) +
    prior(exponential(10), "phi")
  scode <- stancode(bf(y | trials(size)  ~ 1, muy2 ~ x), data = dat,
                         family = dirichlet_multinomial(), prior = prior)
  expect_match2(scode, "array[N, ncat] int Y;")
  expect_match2(scode, "target += dirichlet_multinomial_logit2_lpmf(Y[n] | mu[n], phi);")
  expect_match2(scode, "muy2 += Intercept_muy2 + Xc_muy2 * b_muy2;")
  expect_match2(scode, "lprior += normal_lpdf(b_muy2 | 0, 10);")
  expect_match2(scode, "lprior += cauchy_lpdf(Intercept_muy2 | 0, 1);")
  expect_match2(scode, "lprior += normal_lpdf(Intercept_muy3 | 0, 2);")
  expect_match2(scode, "lprior += exponential_lpdf(phi | 10);")
})

test_that("Stan code for dirichlet models is correct", {
  N <- 15
  dat <- as.data.frame(rdirichlet(N, c(3, 2, 1)))
  names(dat) <- c("y1", "y2", "y3")
  dat$x <- rnorm(N)
  dat$y <- with(dat, cbind(y1, y2, y3))

  # dirichlet in probability-sum(alpha) concentration
  prior <- prior(normal(0, 5), class = "b", dpar = "muy3") +
    prior(exponential(10), "phi")
  scode <- stancode(bf(y ~ 1, muy3 ~ x), data = dat,
                         family = dirichlet(), prior = prior)
  expect_match2(scode, "array[N] vector[ncat] Y;")
  expect_match2(scode, "target += dirichlet_logit_lpdf(Y[n] | mu[n], phi);")
  expect_match2(scode, "muy3 += Intercept_muy3 + Xc_muy3 * b_muy3;")
  expect_match2(scode, "lprior += normal_lpdf(b_muy3 | 0, 5);")
  expect_match2(scode, "lprior += exponential_lpdf(phi | 10);")

  scode <- stancode(bf(y ~ x, phi ~ x), data = dat,
                         family = dirichlet())
  expect_match2(scode, "target += dirichlet_logit_lpdf(Y[n] | mu[n], phi[n]);")
  expect_match2(scode, "phi += Intercept_phi + Xc_phi * b_phi;")
  expect_match2(scode, "phi = exp(phi);")

  # dirichlet2 in alpha parameterization
  prior <- prior(normal(0, 5), class = "b", dpar = "muy3")
  scode <- stancode(bf(y ~ 1, muy3 ~ x), data = dat,
                         family = brmsfamily("dirichlet2"), prior = prior)
  expect_match2(scode, "array[N] vector[ncat] Y;")
  expect_match2(scode, "muy3 = exp(muy3);")
  expect_match2(scode, "target += dirichlet_lpdf(Y[n] | mu[n]);")
  expect_match2(scode, "muy3 += Intercept_muy3 + Xc_muy3 * b_muy3;")
  expect_match2(scode, "mu[n] = transpose([muy1[n], muy2[n], muy3[n]]);")
  expect_match2(scode, "lprior += normal_lpdf(b_muy3 | 0, 5);")
  expect_match2(scode, "lprior += student_t_lpdf(Intercept_muy1 | 3, 0, 2.5);")
})

test_that("Stan code for logistic_normal models is correct", {
  N <- 15
  dat <- as.data.frame(rdirichlet(N, c(3, 2, 1)))
  names(dat) <- c("y1", "y2", "y3")
  dat$x <- rnorm(N)
  dat$y <- with(dat, cbind(y1, y2, y3))

  prior <- prior(normal(0, 5), class = "b", dpar = "muy3") +
    prior(exponential(10), "sigmay1") +
    prior(lkj(3), "lncor")
  scode <- stancode(bf(y ~ x), data = dat,
                         family = logistic_normal(refcat = "y2"),
                         prior = prior)
  expect_match2(scode, "array[N] vector[ncat] Y;")
  expect_match2(scode, "mu[n] = transpose([muy1[n], muy3[n]]);")
  expect_match2(scode, "vector[ncat-1] sigma = transpose([sigmay1, sigmay3]);")
  expect_match2(scode, "target += logistic_normal_cholesky_cor_lpdf(Y[n] | mu[n], sigma, Llncor, 2);")
  expect_match2(scode, "muy3 += Intercept_muy3 + Xc_muy3 * b_muy3;")
  expect_match2(scode, "lprior += normal_lpdf(b_muy3 | 0, 5);")
  expect_match2(scode, "lprior += exponential_lpdf(sigmay1 | 10);")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(Llncor | 3);")

  prior <- prior(normal(0, 5), class = "b", dpar = "muy3") +
    prior(normal(0, 3), class = "b", dpar = "sigmay2")
  scode <- stancode(bf(y ~ 1, muy3 ~ x, sigmay2 ~ x), data = dat,
                         family = logistic_normal(),
                         prior = prior)
  expect_match2(scode, "array[N] vector[ncat] Y;")
  expect_match2(scode, "mu[n] = transpose([muy2[n], muy3[n]]);")
  expect_match2(scode, "sigma[n] = transpose([sigmay2[n], sigmay3]);")
  expect_match2(scode, "target += logistic_normal_cholesky_cor_lpdf(Y[n] | mu[n], sigma[n], Llncor, 1);")
  expect_match2(scode, "muy3 += Intercept_muy3 + Xc_muy3 * b_muy3;")
  expect_match2(scode, "lprior += normal_lpdf(b_muy3 | 0, 5);")
  expect_match2(scode, "lprior += normal_lpdf(b_sigmay2 | 0, 3);")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(Llncor | 1);")
})

test_that("Stan code for ARMA models is correct", {
  dat <- data.frame(y = rep(1:4, 2), x = 1:8, time = 1:8)
  scode <- stancode(y ~ x + ar(time), dat, student())
  expect_match2(scode, "vector[Kar] ar")
  expect_match2(scode, "err[n] = Y[n] - mu[n];")
  expect_match2(scode, "mu[n] += Err[n, 1:Kar] * ar;")

  scode <- stancode(y ~ x + ma(time, q = 2), dat, student())
  expect_match2(scode, "mu[n] += Err[n, 1:Kma] * ma;")

  expect_warning(
    scode <- stancode(mvbind(y, x) ~ 1, dat, gaussian(),
                           autocor = cor_ar()),
    "Argument 'autocor' should be specified within the 'formula' argument"
  )
  expect_match2(scode, "err_y[n] = Y_y[n] - mu_y[n];")

  bform <- bf(y ~ x, sigma ~ x) + acformula(~arma(time, cov = TRUE))
  scode <- stancode(bform, dat, family = student)
  expect_match2(scode, "student_t_time_het_lpdf(Y | nu, mu, sigma, Lcortime")

  bform <- bf(y ~ exp(eta) - 1, eta ~ x, autocor = ~ar(time), nl = TRUE)
  scode <- stancode(bform, dat, family = student,
                         prior = prior(normal(0, 1), nlpar = eta))
  expect_match2(scode, "mu[n] += Err[n, 1:Kar] * ar;")

  # correlations of latent residuals
  scode <- stancode(
    y ~ x + ar(time, cov = TRUE), dat, family = poisson,
    prior = prior(cauchy(0, 10), class = sderr)
  )
  expect_match2(scode, "Lcortime = cholesky_cor_ar1(ar[1], max_nobs_tg);")
  expect_match2(scode,
    "err = scale_time_err(zerr, sderr, Lcortime, nobs_tg, begin_tg, end_tg);"
  )
  expect_match2(scode, "mu += Intercept + Xc * b + err;")
  expect_match2(scode, "lprior += cauchy_lpdf(sderr | 0, 10)")

  scode <- stancode(
    y ~ x + ar(time), dat, family = poisson,
    prior = prior(cauchy(0, 10), class = sderr)
  )
  expect_match2(scode, "vector<lower=-1,upper=1>[Kar] ar")
  expect_match2(scode, "mu[n] += Err[n, 1:Kar] * ar;")
  expect_match2(scode, "err = sderr * zerr;")
  expect_match2(scode, "mu += Intercept + Xc * b + err;")
  expect_match2(scode, "lprior += cauchy_lpdf(sderr | 0, 10)")

  # apply shrinkage priors on sderr
  scode <- stancode(
    y ~ x + ar(time), dat, family = poisson,
    prior = prior(horseshoe(main = TRUE), class = b) +
      prior(horseshoe(), class = sderr)
  )
  expect_match2(scode, "sderr = scales[(1+Kc):(Kc+1)][1];")
})

test_that("Stan code for compound symmetry models is correct", {
  dat <- data.frame(y = rep(1:4, 2), x = 1:8, time = 1:8)
  scode <- stancode(
    y ~ x + cosy(time), dat,
    prior = prior(normal(0, 2), cosy)
  )
  expect_match2(scode, "real<lower=0,upper=1> cosy;")
  expect_match2(scode, "Lcortime = cholesky_cor_cosy(cosy, max_nobs_tg);")
  expect_match2(scode, "lprior += normal_lpdf(cosy | 0, 2)")

  scode <- stancode(bf(y ~ x + cosy(time), sigma ~ x), dat)
  expect_match2(scode, "normal_time_het_lpdf(Y | mu, sigma, Lcortime")

  scode <- stancode(y ~ x + cosy(time), dat, family = poisson)
  expect_match2(scode, "Lcortime = cholesky_cor_cosy(cosy, max_nobs_tg);")
})

test_that("Stan code for UNSTR covariance terms is correct", {
  dat <- data.frame(y = 1:12, x = rnorm(12), tim = c(5:1, 1:5, c(0, 4)),
                    g = c(rep(3:4, 5), rep(2, 2)))

  scode <- stancode(y ~ x + unstr(tim, g), data = dat)
  expect_match2(scode, "normal_time_hom_flex_lpdf(Y | mu, sigma, Lcortime, nobs_tg, begin_tg, end_tg, Jtime_tg);")
  expect_match2(scode, "cortime[choose(k - 1, 2) + j] = Cortime[j, k];")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(Lcortime | 1);")

  scode <- stancode(
    y ~ x + unstr(tim, g), data = dat,
    family = student(), prior = prior(lkj(4), cortime)
  )
  expect_match2(scode, "student_t_time_hom_flex_lpdf(Y | nu, mu, sigma, Lcortime, nobs_tg, begin_tg, end_tg, Jtime_tg);")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(Lcortime | 4);")

  # test standard error
  scode <- stancode(
    y | se(1, sigma = TRUE) ~ x + unstr(tim, g),
    data = dat, family = gaussian(),
  )
  expect_match2(scode, "normal_time_hom_se_flex_lpdf(Y | mu, sigma, se2, Lcortime, nobs_tg, begin_tg, end_tg, Jtime_tg);")

  # test latent representation
  scode <- stancode(
    y ~ x + unstr(tim, g), data = dat,
    family = poisson()
  )
  expect_match2(scode, "err = scale_time_err_flex(zerr, sderr, Lcortime, nobs_tg, begin_tg, end_tg,")
  expect_match2(scode, "mu += Intercept + Xc * b + err;")

  # non-linear model
  scode <- stancode(
    bf(y ~ a, a ~ x, autocor = ~ unstr(tim, g), nl = TRUE),
    data = dat, family = student(), prior = prior(normal(0,1), nlpar = a)
  )
  expect_match2(scode, "student_t_time_hom_flex_lpdf(Y | nu, mu, sigma, Lcortime, nobs_tg, begin_tg, end_tg, Jtime_tg);")
})

test_that("Stan code for intercept only models is correct", {
  expect_match2(stancode(rating ~ 1, data = inhaler),
               "b_Intercept = Intercept;")
  expect_match2(stancode(rating ~ 1, data = inhaler, family = cratio()),
               "b_Intercept = Intercept;")
  expect_match2(stancode(rating ~ 1, data = inhaler, family = categorical()),
               "b_mu3_Intercept = Intercept_mu3;")
})

test_that("Stan code of ordinal models is correct", {
  dat <- data.frame(y = c(rep(1:4, 2), 1, 1), x1 = rnorm(10),
                    x2 = rnorm(10), g = factor(rep(1:2, 5)))

  scode <- stancode(
    y ~ x1, dat, family = cumulative("logit"),
    prior = prior(normal(0, 2), Intercept, coef = 2)
  )
  expect_match2(scode, "target += ordered_logistic_glm_lpmf(Y | Xc, b, Intercept);")
  expect_match2(scode, "lprior += student_t_lpdf(Intercept[1] | 3, 0, 2.5);")
  expect_match2(scode, "lprior += normal_lpdf(Intercept[2] | 0, 2);")

  scode <- stancode(
    y ~ x1, dat, cumulative("probit", threshold = "equidistant"),
    prior = prior(normal(0, 2), Intercept)
  )
  expect_match2(scode, "real cumulative_probit_lpmf(int y")
  expect_match2(scode, "p = Phi(disc * (thres[1] - mu));")
  expect_match2(scode, "real<lower=0> delta;")
  expect_match2(scode, "Intercept[k] = first_Intercept + (k - 1.0) * delta;")
  expect_match2(scode, "b_Intercept = Intercept + dot_product(means_X, b);")
  expect_match2(scode, "lprior += normal_lpdf(first_Intercept | 0, 2);")
  expect_match2(scode, "target += ordered_probit_lpmf(Y[n] | mu[n], Intercept);")

  scode <- stancode(y ~ x1, dat, family = cratio("probit"))
  expect_match2(scode, "real cratio_probit_lpmf(int y")
  expect_match2(scode, "q[k] = std_normal_lcdf(disc * (mu - thres[k]));")

  scode <- stancode(y ~ x1 + cs(x2) + cs(g), dat, family = sratio())
  expect_match2(scode, "real sratio_logit_lpmf(int y")
  expect_match2(scode, "matrix[N, Kcs] Xcs;")
  expect_match2(scode, "matrix[Kcs, nthres] bcs;")
  expect_match2(scode, "mucs = Xcs * bcs;")
  expect_match2(scode,
    "target += sratio_logit_lpmf(Y[n] | mu[n], disc, Intercept - transpose(mucs[n]));"
  )

  scode <- stancode(y ~ x1 + cse(x2) + (cse(1)|g), dat, family = acat())
  expect_match2(scode, "real acat_logit_lpmf(int y")
  expect_match2(scode, "mucs[n, 1] = mucs[n, 1] + r_1_1[J_1[n]] * Z_1_1[n];")
  expect_match2(scode, "b_Intercept = Intercept + dot_product(means_X, b);")

  scode <- stancode(y ~ x1 + (cse(x2)||g), dat, family = acat("probit_approx"))
  expect_match2(scode,
    paste("mucs[n, 3] = mucs[n, 3] + r_1_3[J_1[n]] * Z_1_3[n]",
          "+ r_1_6[J_1[n]] * Z_1_6[n];"))
  expect_match2(scode,
    "target += acat_probit_approx_lpmf(Y[n] | mu[n], disc, Intercept - transpose(mucs[n]));"
  )

  # sum-to-zero thresholds
  scode <- stancode(
    y ~ x1, dat, cumulative("cloglog", threshold = "sum_to_zero"),
    prior = prior(normal(0, 2), Intercept)
  )
  expect_match2(scode, "Intercept_stz = Intercept - mean(Intercept);")
  expect_match2(scode, "cumulative_cloglog_lpmf(Y[n] | mu[n], disc, Intercept_stz);")
  expect_match2(scode, "vector[nthres] b_Intercept = Intercept_stz;")

  # non-linear ordinal models
  scode <- stancode(
    bf(y ~ eta, eta ~ x1, nl = TRUE), dat, family = cumulative(),
    prior = prior(normal(0, 2), nlpar = eta)
  )
  expect_match2(scode, "ordered[nthres] Intercept;")
  expect_match2(scode,
    "target += ordered_logistic_lpmf(Y[n] | mu[n], Intercept);"
  )

  # ordinal mixture models with fixed intercepts
  scode <- stancode(
    bf(y ~ 1, mu1 ~ x1, mu2 ~ 1), data = dat,
    family = mixture(cumulative(), nmix = 2, order = "mu")
  )
  expect_match2(scode, "Intercept_mu2 = fixed_Intercept;")
  expect_match2(scode, "lprior += student_t_lpdf(fixed_Intercept | 3, 0, 2.5);")
})

test_that("Stan code for xbeta models is correct", {
  dat <- data.frame(y = rbeta(10, 1, 1), x = rnorm(10))
  scode <- stancode(y ~ x, dat, family = xbeta())
  expect_match2(scode, "target += xbeta_lpdf(Y | mu, phi, kappa);")
  expect_match2(scode, "lprior += gamma_lpdf(kappa | 0.01, 0.01);")
})

test_that("ordinal disc parameters appear in the Stan code", {
  scode <- stancode(
    bf(rating ~ period + carry + treat, disc ~ period),
    data = inhaler, family = cumulative(),
    prior = prior(normal(0,5), dpar = disc)
  )
  expect_match2(scode,
    "target += cumulative_logit_lpmf(Y[n] | mu[n], disc[n], Intercept)"
  )
  expect_match2(scode, "lprior += normal_lpdf(b_disc | 0, 5)")
  expect_match2(scode, "disc = exp(disc)")
})

test_that("grouped ordinal thresholds appear in the Stan code", {
  dat <- data.frame(
    y = sample(1:6, 10, TRUE),
    y2 = sample(1:6, 10, TRUE),
    gr = rep(c("a", "b"), each = 5),
    th = rep(5:6, each = 5),
    x = rnorm(10)
  )

  prior <- prior(normal(0,1), class = "Intercept", group = "b")
  scode <- stancode(
    y | thres(th, gr) ~ x, data = dat,
    family = sratio(), prior = prior
  )
  expect_match2(scode, "array[ngrthres] int<lower=1> nthres;")
  expect_match2(scode, "merged_Intercept[Kthres_start[1]:Kthres_end[1]] = Intercept_1;")
  expect_match2(scode, "target += sratio_logit_merged_lpmf(Y[n]")
  expect_match2(scode, "lprior += normal_lpdf(Intercept_2 | 0, 1);")
  # centering needs to be deactivated automatically
  expect_match2(scode, "vector[nthres[1]] b_Intercept_1 = Intercept_1;")

  # model with equidistant thresholds
  scode <- stancode(
    y | thres(th, gr) ~ x, data = dat,
    family = cumulative(threshold = "equidistant"),
    prior = prior
  )
  expect_match2(scode, "target += ordered_logistic_merged_lpmf(Y[n]")
  expect_match2(scode, "real first_Intercept_1;")
  expect_match2(scode, "lprior += normal_lpdf(first_Intercept_2 | 0, 1);")
  expect_match2(scode, "Intercept_2[k] = first_Intercept_2 + (k - 1.0) * delta_2;")

  # sum-to-zero constraints
  scode <- stancode(
    y | thres(gr = gr) ~ x, data = dat,
    cumulative(threshold = "sum_to_zero"),
    prior = prior(normal(0, 2), Intercept)
  )
  expect_match2(scode, "merged_Intercept_stz[Kthres_start[2]:Kthres_end[2]] = Intercept_stz_2;")
  expect_match2(scode, "ordered_logistic_merged_lpmf(Y[n] | mu[n], merged_Intercept_stz, Jthres[n]);")

  # ordinal mixture model
  scode <- stancode(
    y | thres(th, gr) ~ x, data = dat,
    family = mixture(cratio, acat, order = "mu"),
    prior = prior
  )
  expect_match2(scode, "ps[1] = log(theta1) + cratio_logit_merged_lpmf(Y[n]")
  expect_match2(scode, "ps[2] = log(theta2) + acat_logit_merged_lpmf(Y[n]")
  expect_match2(scode, "vector[nmthres] merged_Intercept_mu1;")
  expect_match2(scode, "merged_Intercept_mu2[Kthres_start[1]:Kthres_end[1]] = Intercept_mu2_1;")
  expect_match2(scode, "vector[nthres[1]] b_mu1_Intercept_1 = Intercept_mu1_1;")

  # multivariate ordinal model
  bform <- bf(y | thres(th, gr) ~ x, family = sratio) +
    bf(y2 | thres(th, gr) ~ x, family = cumulative)
  scode <- stancode(bform, data = dat)
  expect_match2(scode, "lprior += student_t_lpdf(Intercept_y2_1 | 3, 0, 2.5);")
  expect_match2(scode, "merged_Intercept_y[Kthres_start_y[2]:Kthres_end_y[2]] = Intercept_y_2;")
})


test_that("Stan code of hurdle cumulative model is correct", {
  dat <- data.frame(y = rep(0:4, 2),
                    x1 = rnorm(10),
                    x2 = rnorm(10),
                    g = factor(rep(1:2, 5)))

  scode <- stancode(
    y ~ x1, dat, family = hurdle_cumulative("logit"),
    prior = prior(normal(0, 2), Intercept, coef = 2)
  )
  expect_match2(scode,
    "target += hurdle_cumulative_ordered_logistic_lpmf(Y[n] | mu[n], hu, disc, Intercept);"
  )

  scode <- stancode(
    bf(y ~ x1, hu ~ x2), dat,
    hurdle_cumulative("probit", threshold = "equidistant"),
    prior = prior(normal(0, 2), Intercept)
  )

  expect_match2(scode, "real hurdle_cumulative_ordered_probit_lpmf(int y")
  expect_match2(scode, "p = Phi(disc * (thres[1] - mu)) * (1 - hu);")
  expect_match2(scode, "Intercept[k] = first_Intercept + (k - 1.0) * delta;")

  # sum-to-zero thresholds
  scode <- stancode(
    bf(y ~ x1, hu ~ x2, disc ~ g), dat,
    hurdle_cumulative("cloglog", threshold = "sum_to_zero"),
    prior = prior(normal(0, 2), Intercept)
  )
  expect_match2(scode, "Intercept_stz = Intercept - mean(Intercept);")
  expect_match2(scode, "hurdle_cumulative_cloglog_lpmf(Y[n] | mu[n], hu[n], disc[n], Intercept_stz);")
  expect_match2(scode, "vector[nthres] b_Intercept = Intercept_stz;")


  # non-linear ordinal models
  scode <- stancode(
    bf(y ~ eta, eta ~ x1, nl = TRUE),
    dat,
    family = hurdle_cumulative(),
    prior = prior(normal(0, 2), nlpar = eta)
  )
  expect_match2(scode,
    "target += hurdle_cumulative_ordered_logistic_lpmf(Y[n] | mu[n], hu, disc, Intercept);"
  )
})

test_that("monotonic effects appear in the Stan code", {
  dat <- data.frame(y = rpois(120, 10), x1 = rep(1:4, 30),
                    x2 = factor(rep(c("a", "b", "c"), 40), ordered = TRUE),
                    g = rep(1:10, each = 12))

  prior <- c(prior(normal(0,1), class = b, coef = mox1),
             prior(dirichlet(c(1,0.5,2)), simo, coef = mox11),
             prior(dirichlet(c(1,0.5,2)), simo, coef = mox21))
  scode <- stancode(y ~ y*mo(x1)*mo(x2), dat, prior = prior)
  expect_match2(scode, "array[N] int Xmo_3;")
  expect_match2(scode, "simplex[Jmo[1]] simo_1;")
  expect_match2(scode, "(bsp[2]) * mo(simo_2, Xmo_2[n])")
  expect_match2(scode,
    "(bsp[6]) * mo(simo_7, Xmo_7[n]) * mo(simo_8, Xmo_8[n]) * Csp_3[n]"
  )
  expect_match2(scode, "lprior += normal_lpdf(bsp[1] | 0, 1)")
  expect_match2(scode, "lprior += dirichlet_lpdf(simo_1 | con_simo_1);")
  expect_match2(scode, "lprior += dirichlet_lpdf(simo_8 | con_simo_8);")

  scode <- stancode(y ~ mo(x1) + (mo(x1) | x2), dat)
  expect_match2(scode, "(bsp[1] + r_1_2[J_1[n]]) * mo(simo_1, Xmo_1[n])")
  expect_true(!grepl("Z_1_w", scode))

  # test issue reported in discourse post #12978
  scode <- stancode(y ~ mo(x1) + (mo(x1) | x2) + (mo(x1) | g), dat)
  expect_match2(scode, "(bsp[1] + r_1_2[J_1[n]] + r_2_2[J_2[n]]) * mo(simo_1, Xmo_1[n])")

  # test issue #813
  scode <- stancode(y ~ mo(x1):y, dat)
  expect_match2(scode, "mu[n] += (bsp[1]) * mo(simo_1, Xmo_1[n]) * Csp_1[n];")

  # test issue #924 (conditional monotonicity)
  prior <- c(prior(dirichlet(c(1,0.5,2)), simo, coef = "v"),
             prior(dirichlet(c(1,0.5,2)), simo, coef = "w"))
  scode <- stancode(y ~ y*mo(x1, id = "v")*mo(x2, id = "w"),
                         dat, prior = prior)
  expect_match2(scode, "lprior += dirichlet_lpdf(simo_1 | con_simo_1);")
  expect_match2(scode, "lprior += dirichlet_lpdf(simo_2 | con_simo_2);")
  expect_match2(scode, "simplex[Jmo[6]] simo_6 = simo_2;")
  expect_match2(scode, "simplex[Jmo[7]] simo_7 = simo_1;")

  expect_error(
    stancode(y ~ mo(x1) + (mo(x2) | x2), dat),
    "Special group-level terms require"
  )

  prior <- prior(beta(1, 1), simo, coef = mox11)
  expect_error(
    stancode(y ~ mo(x1), dat, prior = prior),
    "'dirichlet' is the only valid prior for simplex parameters"
  )
})

test_that("Stan code for non-linear models is correct", {
  flist <- list(a ~ x, b ~ z + (1|g))
  data <- data.frame(
    y = rgamma(9, 1, 1), x = rnorm(9),
    z = rnorm(9), v = 1L:9L, g = rep(1:3, 3)
  )
  prior <- c(set_prior("normal(0,5)", nlpar = "a"),
             set_prior("normal(0,1)", nlpar = "b"))
  # syntactic validity is already checked within stancode
  scode <- stancode(
    bf(y ~ a - exp(b^z) * (z <= a) * v, flist = flist, nl = TRUE),
    data = data, prior = prior
  )
  expect_match2(scode,
    "mu[n] = (nlp_a[n] - exp(nlp_b[n] ^ C_1[n]) * (C_1[n] <= nlp_a[n]) * C_2[n]);"
  )
  expect_match2(scode, "vector[N] C_1;")
  expect_match2(scode, "array[N] int C_2;")

  # non-linear predictor can be computed outside a loop
  scode <- stancode(bf(y ~ a - exp(b + z), flist = flist,
                            nl = TRUE, loop = FALSE),
                         data = data, prior = prior)
  expect_match2(scode, "mu = (nlp_a - exp(nlp_b + C_1));")

  # check if that also works with threading
  scode <- stancode(bf(y ~ a - exp(b + z), flist = flist,
                            nl = TRUE, loop = FALSE),
                         data = data, prior = prior,
                         threads = threading(2), parse = FALSE)
  expect_match2(scode, "mu = (nlp_a - exp(nlp_b + C_1[start:end]));")


  flist <- list(a1 ~ 1, a2 ~ z + (x|g))
  prior <- c(set_prior("beta(1,1)", nlpar = "a1", lb = 0, ub = 1),
             set_prior("normal(0,1)", nlpar = "a2"))
  scode <- stancode(
    bf(y ~ a1 * exp(-x/(a2 + z)),
       flist = flist, nl = TRUE),
    data = data, family = Gamma("log"),
    prior = prior
  )
  expect_match2(scode, "mu[n] = exp(nlp_a1[n] * exp( - C_1[n] / (nlp_a2[n] + C_2[n])));")

  bform <- bf(y ~ x) +
    nlf(sigma ~ a1 * exp(-x/(a2 + z))) +
    lf(a1 ~ 1, a2 ~ z + (x|g)) +
    lf(alpha ~ x)
  scode <- stancode(
    bform, data, family = skew_normal(),
    prior = c(
      prior(normal(0, 1), nlpar = a1),
      prior(normal(0, 5), nlpar = a2)
    )
  )
  expect_match2(scode, "nlp_a1 += X_a1 * b_a1")
  expect_match2(scode,
    "sigma[n] = exp(nlp_a1[n] * exp( - C_sigma_1[n] / (nlp_a2[n] + C_sigma_2[n])))"
  )
  expect_match2(scode, "lprior += normal_lpdf(b_a2 | 0, 5)")
})

test_that("Stan code for nested non-linear parameters is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), z = 1:5)
  bform <- bf(
    y ~ lb + (1 - lb) * inv_logit(b * x),
    b + a ~ 1 + (1 | z), nlf(lb ~ inv_logit(a / x)),
    nl = TRUE
  )
  bprior <- prior(normal(0, 1), nlpar = "a") +
    prior(normal(0, 1), nlpar = "b")
  scode <- stancode(bform, dat, prior = bprior)
  expect_match2(scode, "nlp_lb[n] = (inv_logit(nlp_a[n] / C_lb_1[n]));")
  expect_match2(scode,
    "mu[n] = (nlp_lb[n] + (1 - nlp_lb[n]) * inv_logit(nlp_b[n] * C_1[n]));"
  )
})

test_that("stancode is correct for non-linear matrix covariates", {
  N <- 10
  dat <- data.frame(y=rnorm(N))
  dat$X <- matrix(rnorm(N*2), N, 2)
  dat$X2 <- matrix(1L:4L, N, 2)

  # numeric matrix
  nlfun_stan <- "
    real nlfun(real a, real b, real c, row_vector X) {
       return a + b * X[1] + c * X[2];
    }
  "
  nlstanvar <- stanvar(scode = nlfun_stan, block = "functions")
  bform <- bf(y~nlfun(a, b, c, X), a~1, b~1, c~1, nl = TRUE)
  scode <- stancode(bform, dat, stanvars = nlstanvar)
  expect_match2(scode, "matrix[N, 2] C_1;")

  # integer matrix
  nlfun_stan_int <- "
    real nlfun(real a, real b, real c, array[] int X) {
       return a + b * X[1] + c * X[2];
    }
  "
  nlstanvar <- stanvar(scode = nlfun_stan_int, block = "functions")
  bform <- bf(y~nlfun(a, b, c, X2), a~1, b~1, c~1, nl = TRUE)
  scode <- stancode(bform, dat, stanvars = nlstanvar)
  expect_match2(scode, "array[N, 2] int C_1;")
})

test_that("stancode accepts very long non-linear formulas", {
  data <- data.frame(y = rnorm(10), this_is_a_very_long_predictor = rnorm(10))
  expect_silent(stancode(bf(y ~ b0 + this_is_a_very_long_predictor +
                                 this_is_a_very_long_predictor +
                                 this_is_a_very_long_predictor,
                                 b0 ~ 1, nl = TRUE),
                data = data, prior = prior(normal(0,1), nlpar = "b0")))
})

test_that("no loop in trans-par is defined for simple 'identity' models", {
  expect_true(!grepl(stancode(time ~ age, data = kidney),
                     "mu[n] = (mu[n]);", fixed = TRUE))
  expect_true(!grepl(stancode(time ~ age, data = kidney,
                                   family = poisson("identity")),
                     "mu[n] = (mu[n]);", fixed = TRUE))
})

test_that("known standard errors appear in the Stan code", {
  scode <- stancode(time | se(age) ~ sex, data = kidney)
  expect_match2(scode, "target += normal_lpdf(Y | mu, se)")
  scode <- stancode(time | se(age) + weights(age) ~ sex, data = kidney)
  expect_match2(scode, "target += weights[n] * (normal_lpdf(Y[n] | mu[n], se[n]))")
  scode <- stancode(time | se(age, sigma = TRUE) ~ sex, data = kidney)
  expect_match2(scode, "target += normal_lpdf(Y | mu, sqrt(square(sigma) + se2))")
  scode <- stancode(bf(time | se(age, sigma = TRUE) ~ sex, sigma ~ sex),
                         data = kidney)
  expect_match2(scode, "target += normal_lpdf(Y | mu, sqrt(square(sigma) + se2))")
})

test_that("functions defined in 'stan_funs' appear in the functions block", {
  test_fun <- paste0("  real test_fun(real a, real b) {\n",
                     "    return a + b;\n",
                     "  }\n")
  scode <- SW(stancode(time ~ age, data = kidney, stan_funs = test_fun))
  expect_match2(scode, test_fun)
})

test_that("FCOR matrices appear in the Stan code", {
  data <- data.frame(y = 1:5)
  V <- diag(5)
  expect_match2(stancode(y ~ fcor(V), data = data, family = gaussian(),
                              data2 = list(V = V)),
               "target += normal_fcor_hom_lpdf(Y | mu, sigma, Lfcor);")
  expect_match2(stancode(y ~ fcor(V), data = data, family = student(),
                              data2 = list(V = V)),
               "target += student_t_fcor_hom_lpdf(Y | nu, mu, sigma, Lfcor);")
})

test_that("Stan code for GAMMs is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = factor(rep(1:2, 5)))
  scode <- stancode(y ~ s(x) + (1|g), data = dat,
                         prior = set_prior("normal(0,2)", "sds"))
  expect_match2(scode, "Zs_1_1 * s_1_1")
  expect_match2(scode, "matrix[N, knots_1[1]] Zs_1_1")
  expect_match2(scode, "target += std_normal_lpdf(zs_1_1)")
  expect_match2(scode, "lprior += normal_lpdf(sds_1 | 0,2)")

  prior <- c(set_prior("normal(0,5)", nlpar = "lp"),
             set_prior("normal(0,2)", "sds", nlpar = "lp"))
  scode <- stancode(bf(y ~ lp, lp ~ s(x) + (1|g), nl = TRUE),
                         data = dat, prior = prior)
  expect_match2(scode, "Zs_lp_1_1 * s_lp_1_1")
  expect_match2(scode, "matrix[N, knots_lp_1[1]] Zs_lp_1_1")
  expect_match2(scode, "target += std_normal_lpdf(zs_lp_1_1)")
  expect_match2(scode, "lprior += normal_lpdf(sds_lp_1 | 0,2)")

  scode <- stancode(
    y ~ s(x) + t2(x,y), data = dat,
    prior = set_prior("normal(0,1)", "sds") +
      set_prior("normal(0,2)", "sds", coef = "t2(x, y)")
  )
  expect_match2(scode, "Zs_2_2 * s_2_2")
  expect_match2(scode, "matrix[N, knots_2[2]] Zs_2_2")
  expect_match2(scode, "target += std_normal_lpdf(zs_2_2)")
  expect_match2(scode, "lprior += normal_lpdf(sds_1 | 0,1)")
  expect_match2(scode, "lprior += normal_lpdf(sds_2 | 0,2)")

  scode <- stancode(y ~ g + s(x, by = g), data = dat)
  expect_match2(scode, "vector[knots_2[1]] zs_2_1")
  expect_match2(scode, "s_2_1 = sds_2[1] * zs_2_1")
})

test_that("Stan code of response times models is correct", {
  dat <- epilepsy
  dat$cens <- sample(-1:1, nrow(dat), TRUE)
  scode <- stancode(count ~ Trt + (1|patient),
                      data = dat, family = exgaussian("log"),
                      prior = prior(gamma(1,1), class = beta))
  expect_match2(scode,
    "target += exp_mod_normal_lpdf(Y | mu - beta, sigma, inv(beta))"
  )
  expect_match2(scode, "mu = exp(mu)")
  expect_match2(scode, "lprior += gamma_lpdf(beta | 1, 1)")

  scode <- stancode(bf(count ~ Trt + (1|patient),
                         sigma ~ Trt, beta ~ Trt),
                      data = dat, family = exgaussian())
  expect_match2(scode,
    "target += exp_mod_normal_lpdf(Y | mu - beta, sigma, inv(beta))"
  )
  expect_match2(scode, "beta = exp(beta)")

  scode <- stancode(count | cens(cens) ~ Trt + (1|patient),
                      data = dat, family = exgaussian("inverse"))
  expect_match2(scode,
    "target += exp_mod_normal_lccdf(Y[Jrcens[1:Nrcens]] | mu[Jrcens[1:Nrcens]] - beta, sigma, inv(beta));"
  )

  scode <- stancode(count ~ Trt, dat, family = shifted_lognormal())
  expect_match2(scode, "target += lognormal_lpdf(Y - ndt | mu, sigma)")

  scode <- stancode(count | cens(cens) ~ Trt, dat, family = shifted_lognormal())
  expect_match2(scode,
    "target += lognormal_lcdf(Y[Jlcens[1:Nlcens]] - ndt | mu[Jlcens[1:Nlcens]], sigma);"
  )

  # test issue #837
  scode <- stancode(mvbind(count, zBase) ~ Trt, data = dat,
                         family = shifted_lognormal())
  expect_match2(scode, "lprior += uniform_lpdf(ndt_count | 0, min_Y_count)")
  expect_match2(scode, "lprior += uniform_lpdf(ndt_zBase | 0, min_Y_zBase)")
})

test_that("Stan code of wiener diffusion models is correct", {
  dat <- data.frame(q = 1:10, resp = sample(0:1, 10, TRUE), x = rnorm(10))
  scode <- stancode(q | dec(resp) ~ x, data = dat, family = wiener())
  expect_match2(scode,
    "target += wiener_diffusion_lpdf(Y[n] | dec[n], bs, ndt, bias, mu[n])"
  )

  scode <- stancode(bf(q | dec(resp) ~ x, bs ~ x, ndt ~ x, bias ~ x),
                         data = dat, family = wiener())
  expect_match2(scode,
    "target += wiener_diffusion_lpdf(Y[n] | dec[n], bs[n], ndt[n], bias[n], mu[n])"
  )
  expect_match2(scode, "bias = inv_logit(bias);")

  scode <- stancode(bf(q | dec(resp) ~ x, ndt = 0.5),
                         data = dat, family = wiener())
  expect_match2(scode, "real ndt = 0.5;")

  expect_error(stancode(q ~ x, data = dat, family = wiener()),
               "Addition argument 'dec' is required for family 'wiener'")
})

test_that("Group IDs appear in the Stan code", {
  form <- bf(count ~ Trt + (1+Trt|3|visit) + (1|patient),
             shape ~ (1|3|visit) + (Trt||patient))
  scode <- stancode(form, data = epilepsy, family = negbinomial())
  expect_match2(scode, "r_2_1 = r_2[, 1]")
  expect_match2(scode, "r_2_shape_3 = r_2[, 3]")

  form <- bf(count ~ a, sigma ~ (1|3|visit) + (Trt||patient),
             a ~ Trt + (1+Trt|3|visit) + (1|patient), nl = TRUE)
  scode <- stancode(form, data = epilepsy, family = student(),
                      prior = set_prior("normal(0,5)", nlpar = "a"))
  expect_match2(scode, "r_2_a_2 = r_2[, 2];")
  expect_match2(scode, "r_1_sigma_2 = (sd_1[2] * (z_1[2]));")
})

test_that("weighted, censored, and truncated likelihoods are correct", {
  dat <- data.frame(y = 1:9, x = rep(-1:1, 3), y2 = 10:18)

  scode <- stancode(y | weights(y2) ~ 1, dat, poisson())
  expect_match2(scode, "target += weights[n] * (poisson_log_lpmf(Y[n] | mu[n]));")

  scode <- stancode(y | trials(y2) + weights(y2) ~ 1, dat, binomial())
  expect_match2(scode,
    "target += weights[n] * (binomial_logit_lpmf(Y[n] | trials[n], mu[n]));"
  )

  scode <- stancode(y | cens(x, y2) ~ 1, dat, family = poisson())
  expect_match2(scode, "target += poisson_lpmf(Y[n] | mu[n]);")
  expect_match2(scode, "poisson_lcdf(rcens[n] | mu[n])")

  scode <- stancode(y | cens(x) ~ 1, dat, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lccdf(Y[n] | mu[n], sigma, quantile);")

  scode <- stancode(bf(y | cens(x) ~ 1, shape ~ 1), dat, family = Gamma())
  expect_match2(scode, "target += gamma_lpdf(Y[Jevent[1:Nevent]] | shape[Jevent[1:Nevent]], shape[Jevent[1:Nevent]] ./ mu[Jevent[1:Nevent]]);")

  dat$x[1] <- 2
  scode <- stancode(y | cens(x, y2) ~ 1, dat, family = asym_laplace())
  expect_match2(scode, "target += log_diff_exp(\n")
  expect_match2(scode, "asym_laplace_lcdf(rcens[n] | mu[n], sigma, quantile),")

  dat$x <- 1
  expect_match2(stancode(y | cens(x) + weights(x) ~ 1, dat, exponential()),
   "target += weights[n] * exponential_lccdf(Y[n] | inv(mu[n]));")

  scode <- stancode(y | cens(x) + trunc(0.1) ~ 1, dat, exponential())
  expect_match2(scode, "target += exponential_lccdf(Y[n] | inv(mu[n])) -")
  expect_match2(scode, "  exponential_lccdf(lb[n] | inv(mu[n]));")

  scode <- stancode(y | cens(x) + trunc(ub = 30) ~ 1, dat)
  expect_match2(scode, "target += normal_lccdf(Y[n] | mu[n], sigma) -")
  expect_match2(scode, "  normal_lcdf(ub[n] | mu[n], sigma);")

  scode <- stancode(y | weights(x) + trunc(0, 30) ~ 1, dat)
  expect_match2(scode, "target += weights[n] * (normal_lpdf(Y[n] | mu[n], sigma) -")
  expect_match2(scode, "  log_diff_exp(normal_lcdf(ub[n] | mu[n], sigma),")

  expect_match2(
    stancode(y | trials(y2) + weights(y2) ~ 1, dat, beta_binomial()),
    "target += weights[n] * (beta_binomial_lpmf(Y[n] | trials[n], mu[n] .* phi,"
  )
  expect_match2(
    stancode(y | trials(y2) + trunc(0, 30) ~ 1, dat, beta_binomial()),
    "log_diff_exp(beta_binomial_lcdf(ub[n] | trials[n], mu[n] .* phi,"
  )
  expect_match2(
    stancode(y | trials(y2) + cens(x, y2) ~ 1, dat, beta_binomial()),
    "beta_binomial_lcdf(rcens[n] | trials[n], mu[n] .* phi,"
  )
})

test_that("noise-free terms appear in the Stan code", {
  N <- 30
  dat <- data.frame(
    y = rnorm(N), x = rnorm(N), z = rnorm(N),
    xsd = abs(rnorm(N, 1)), zsd = abs(rnorm(N, 1)),
    ID = rep(1:5, each = N / 5)
  )

  me_prior <- prior(normal(0,5)) +
    prior(normal(0, 10), "meanme") +
    prior(cauchy(0, 5), "sdme", coef = "mez") +
    prior(lkj(2), "corme")
  scode <- stancode(
    y ~ me(x, xsd)*me(z, zsd)*x, data = dat, prior = me_prior,
    sample_prior = "yes"
  )
  expect_match2(scode,
    "(bsp[1]) * Xme_1[n] + (bsp[2]) * Xme_2[n] + (bsp[3]) * Xme_1[n] * Xme_2[n]"
  )
  expect_match2(scode, "(bsp[6]) * Xme_1[n] * Xme_2[n] * Csp_3[n]")
  expect_match2(scode, "target += normal_lpdf(Xn_2 | Xme_2, noise_2)")
  expect_match2(scode, "lprior += normal_lpdf(bsp | 0, 5)")
  expect_match2(scode, "target += std_normal_lpdf(to_vector(zme_1))")
  expect_match2(scode, "lprior += normal_lpdf(meanme_1 | 0, 10)")
  expect_match2(scode, "lprior += cauchy_lpdf(sdme_1[2] | 0, 5)")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(Lme_1 | 2)")
  expect_match2(scode, "+ transpose(diag_pre_multiply(sdme_1, Lme_1) * zme_1)")
  expect_match2(scode, "corme_1[choose(k - 1, 2) + j] = Corme_1[j, k];")

  scode <- stancode(
    y ~ me(x, xsd)*z + (me(x, xsd)*z | ID), data = dat
  )
  expect_match2(scode, "(bsp[1] + r_1_3[J_1[n]]) * Xme_1[n]")
  expect_match2(scode, "(bsp[2] + r_1_4[J_1[n]]) * Xme_1[n] * Csp_1[n]")

  expect_match2(stancode(y ~ I(me(x, xsd)^2), data = dat),
               "(bsp[1]) * (Xme_1[n]^2)")

  # test that noise-free variables are unique across model parts
  scode <- stancode(
    bf(y ~ me(x, xsd)*me(z, zsd)*x, sigma ~ me(x, xsd)),
    data = dat, prior = prior(normal(0,5))
  )
  expect_match2(scode, "mu[n] += (bsp[1]) * Xme_1[n]")
  expect_match2(scode, "sigma[n] += (bsp_sigma[1]) * Xme_1[n]")

  scode <- stancode(
    bf(y ~ a * b, a + b ~ me(x, xsd), nl = TRUE),
    data = dat,
    prior = prior(normal(0,5), nlpar = a) +
      prior(normal(0, 5), nlpar = b)
  )
  expect_match2(scode, "nlp_a[n] += (bsp_a[1]) * Xme_1[n]")
  expect_match2(scode, "nlp_b[n] += (bsp_b[1]) * Xme_1[n]")

  bform <- bf(mvbind(y, z) ~ me(x, xsd)) +
    set_rescor(TRUE) + set_mecor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode, "mu_y[n] += (bsp_y[1]) * Xme_1[n]")
  expect_match2(scode, "mu_z[n] += (bsp_z[1]) * Xme_1[n]")
  expect_match2(scode, "Xme_1 = meanme_1[1] + sdme_1[1] * zme_1;")

  # noise-free terms with grouping factors
  bform <- bf(y ~ me(x, xsd, ID) + me(z, xsd) + (me(x, xsd, ID) | ID))
  scode <- stancode(bform, dat)
  expect_match2(scode, "vector[Nme_1] Xn_1;")
  expect_match2(scode, "Xme_1 = meanme_1[1] + sdme_1[1] * zme_1;")
  expect_match2(scode, "Xme_2 = meanme_2[1] + sdme_2[1] * zme_2;")
  expect_match2(scode, "(bsp[1] + r_1_2[J_1[n]]) * Xme_1[Jme_1[n]]")

  bform <- bform + set_mecor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode, "Xme_1 = meanme_1[1] + sdme_1[1] * zme_1;")
})

test_that("Stan code of multi-membership models is correct", {
  dat <- data.frame(y = rnorm(10), g1 = sample(1:10, 10, TRUE),
                    g2 = sample(1:10, 10, TRUE), w1 = rep(1, 10),
                    w2 = rep(abs(rnorm(10))))
  expect_match2(stancode(y ~ (1|mm(g1, g2)), data = dat),
    paste0(" W_1_1[n] * r_1_1[J_1_1[n]] * Z_1_1_1[n]",
           " + W_1_2[n] * r_1_1[J_1_2[n]] * Z_1_1_2[n]")
  )
  expect_match2(stancode(y ~ (1+w1|mm(g1,g2)), data = dat),
    paste0(" W_1_1[n] * r_1_2[J_1_1[n]] * Z_1_2_1[n]",
           " + W_1_2[n] * r_1_2[J_1_2[n]] * Z_1_2_2[n]")
  )
  expect_match2(stancode(y ~ (1+mmc(w1, w2)|mm(g1,g2)), data = dat),
    " W_1_2[n] * r_1_2[J_1_2[n]] * Z_1_2_2[n];"
  )
})

test_that("by variables in grouping terms are handled correctly", {
  dat <- data.frame(
    y = rnorm(100), x = rnorm(100),
    g = rep(1:10, each = 10),
    z = factor(rep(c(0, 4.5, 3, 2, 5), each = 20))
  )
  scode <- stancode(y ~ x + (1 | gr(g, by = z)), dat)
  expect_match2(scode, "r_1_1 = (transpose(sd_1[1, Jby_1]) .* (z_1[1]));")
  scode <- stancode(y ~ x + (x | gr(g, by = z)), dat)
  expect_match2(scode, "r_1 = scale_r_cor_by(z_1, sd_1, L_1, Jby_1);")
  expect_match2(scode, "lprior += student_t_lpdf(to_vector(sd_1) | 3, 0, 2.5)")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(L_1[5] | 1);")
})

test_that("Group syntax | and || is handled correctly,", {
  data <- data.frame(y = rnorm(10), x = rnorm(10),
                     g1 = rep(1:5, each = 2), g2 = rep(1:2, 5))
  scode <- stancode(y ~ x + (1+x||g1) + (I(x/4)|g2), data)
  expect_match2(scode, "r_1_2 = (sd_1[2] * (z_1[2]));")
  expect_match2(scode, "r_2_1 = r_2[, 1];")
  expect_match2(scode, "r_2 = scale_r_cor(z_2, sd_2, L_2);")
})

test_that("predicting zi and hu works correctly", {
  scode <- stancode(bf(count ~ Trt, zi ~ Trt), epilepsy,
                         family = "zero_inflated_poisson")
  expect_match2(scode,
    "target += zero_inflated_poisson_log_logit_lpmf(Y[n] | mu[n], zi[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))
  expect_true(!grepl("exp(mu[n])", scode, fixed = TRUE))

  scode <- stancode(bf(count ~ Trt, zi ~ Trt), epilepsy,
                         family = zero_inflated_poisson(identity))
  expect_match2(scode,
    "target += zero_inflated_poisson_logit_lpmf(Y[n] | mu[n], zi[n])"
  )

  scode <- stancode(bf(count ~ Trt, zi ~ Trt), epilepsy,
                         family = "zero_inflated_binomial")
  expect_match2(scode,
    "target += zero_inflated_binomial_blogit_logit_lpmf(Y[n] | trials[n], mu[n], zi[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))

  fam <- zero_inflated_binomial("probit", link_zi = "identity")
  scode <- stancode(
    bf(count ~ Trt, zi ~ Trt), epilepsy, family = fam,
    prior = prior("", class = Intercept, dpar = zi, lb = 0, ub = 1)
  )
  expect_match2(scode,
    "target += zero_inflated_binomial_lpmf(Y[n] | trials[n], mu[n], zi[n])"
  )
  expect_match2(scode, "mu = Phi(mu);")

  scode <- stancode(bf(count ~ Trt, zi ~ Trt), epilepsy,
                         family = "zero_inflated_beta_binomial")
  expect_match2(scode,
                paste("target += zero_inflated_beta_binomial_logit_lpmf(Y[n]",
                      "| trials[n], mu[n], phi, zi[n])"))
  expect_match2(scode, "mu = inv_logit(mu);")
  scode <- stancode(
    bf(count ~ Trt, zi ~ Trt), epilepsy,
    zero_inflated_beta_binomial("probit", link_zi = "identity"),
    prior = prior("", class = Intercept, dpar = zi, lb = 0, ub = 1)
  )
  expect_match2(scode,
                paste("target += zero_inflated_beta_binomial_lpmf(Y[n]",
                      "| trials[n], mu[n], phi, zi[n])"))
  expect_match2(scode, "mu = Phi(mu);")

  scode <- stancode(
    bf(count ~ Trt, zi ~ Trt), epilepsy,
    family = zero_inflated_beta()
  )
  expect_match2(scode,
    "target += zero_inflated_beta_logit_lpdf(Y[n] | mu[n], phi, zi[n])"
  )

  scode <- stancode(bf(count ~ Trt, hu ~ Trt), epilepsy,
                         family = "hurdle_negbinomial")
  expect_match2(scode,
    "target += hurdle_neg_binomial_log_logit_lpmf(Y[n] | mu[n], shape, hu[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))
  expect_true(!grepl("exp(mu)", scode, fixed = TRUE))

  scode <- stancode(bf(count ~ Trt, hu ~ Trt), epilepsy,
                         family = "hurdle_gamma")
  expect_match2(scode,
    "hurdle_gamma_logit_lpdf(Y[n] | shape, shape / mu[n], hu[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))

  scode <- stancode(
    bf(count ~ Trt, hu ~ Trt), epilepsy,
    family = hurdle_gamma(link_hu = "identity"),
    prior = prior("", class = Intercept, dpar = hu, lb = 0, ub = 1)
  )
  expect_match2(scode, "hurdle_gamma_lpdf(Y[n] | shape, shape / mu[n], hu[n])")
  expect_true(!grepl("inv_logit\\(", scode))
})

test_that("fixing auxiliary parameters is possible", {
  scode <- stancode(bf(y ~ 1, sigma = 0.5), data = list(y = rnorm(10)))
  expect_match2(scode, "real sigma = 0.5;")
})

test_that("Stan code of quantile regression models is correct", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  scode <- stancode(y ~ x, data, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lpdf(Y[n] | mu[n], sigma, quantile)")

  scode <- stancode(bf(y ~ x, quantile = 0.75), data, family = asym_laplace())
  expect_match2(scode, "real quantile = 0.75;")

  scode <- stancode(y | cens(c) ~ x, data, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lccdf(Y[n] | mu[n], sigma, quantile)")

  scode <- stancode(bf(y ~ x, sigma ~ x), data, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lpdf(Y[n] | mu[n], sigma[n], quantile)")

  scode <- stancode(bf(y ~ x, quantile = 0.75), data,
                         family = brmsfamily("zero_inflated_asym_laplace"))
  expect_match2(scode,
    "target += zero_inflated_asym_laplace_lpdf(Y[n] | mu[n], sigma, quantile, zi)"
  )
})

test_that("Stan code of addition term 'rate' is correct", {
  data <- data.frame(y = rpois(10, 1), x = rnorm(10), time = 1:10)
  scode <- stancode(y | rate(time) ~ x, data, poisson())
  expect_match2(scode, "target += poisson_log_lpmf(Y | mu + log_denom);")

  scode <- stancode(y | rate(time) ~ x, data, poisson("identity"))
  expect_match2(scode, "target += poisson_lpmf(Y | mu .* denom);")

  scode <- stancode(y | rate(time) ~ x, data, negbinomial())
  expect_match2(scode, "target += neg_binomial_2_log_lpmf(Y | mu + log_denom, shape .* denom);")

  bform <- bf(y | rate(time) ~ mi(x), shape ~ mi(x), family = negbinomial()) +
    bf(x | mi() ~ 1, family = gaussian())
  scode <- stancode(bform, data)
  expect_match2(scode, "target += neg_binomial_2_log_lpmf(Y_y | mu_y + log_denom_y, shape_y .* denom_y);")

  scode <- stancode(y | rate(time) ~ x, data, brmsfamily("negbinomial2"))
  expect_match2(scode, "target += neg_binomial_2_log_lpmf(Y | mu + log_denom, inv(sigma) .* denom);")

  scode <- stancode(y | rate(time) + cens(1) ~ x, data, geometric())
  expect_match2(scode,
    "target += neg_binomial_2_lpmf(Y[Jevent[1:Nevent]] | mu[Jevent[1:Nevent]] .* denom[Jevent[1:Nevent]], 1 .* denom[Jevent[1:Nevent]]);"
  )
})

test_that("Stan code of GEV models is correct", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  SW(scode <- stancode(y ~ x, data, gen_extreme_value()))
  expect_match2(scode, "target += gen_extreme_value_lpdf(Y[n] | mu[n], sigma, xi)")
  expect_match2(scode, "xi = scale_xi(tmp_xi, Y, mu, sigma)")

  SW(scode <- stancode(bf(y ~ x, sigma ~ x), data, gen_extreme_value()))
  expect_match2(scode, "xi = scale_xi(tmp_xi, Y, mu, sigma)")

  SW(scode <- stancode(bf(y ~ x, xi ~ x), data, gen_extreme_value()))
  expect_match2(scode, "xi = expm1(xi)")

  SW(scode <- stancode(bf(y ~ x, xi = 0), data, gen_extreme_value()))
  expect_match2(scode, "real xi = 0;  // shape parameter")

  SW(scode <- stancode(y | cens(c) ~ x, data, gen_extreme_value()))
  expect_match2(scode, "target += gen_extreme_value_lccdf(Y[n] | mu[n], sigma, xi)")
})

test_that("Stan code of Cox models is correct", {
  data <- data.frame(y = rexp(100), ce = sample(0:1, 100, TRUE),
                     x = rnorm(100), g = sample(1:3, 100, TRUE))
  bform <- bf(y | cens(ce) ~ x)
  scode <- stancode(bform, data, brmsfamily("cox"))
  expect_match2(scode,
    "target += cox_log_lpdf(Y[Jevent[1:Nevent]] | mu[Jevent[1:Nevent]], bhaz[Jevent[1:Nevent]], cbhaz[Jevent[1:Nevent]]);"
  )
  expect_match2(scode, "vector[N] cbhaz = Zcbhaz * sbhaz;")
  expect_match2(scode, "lprior += dirichlet_lpdf(sbhaz | con_sbhaz);")
  expect_match2(scode, "simplex[Kbhaz] sbhaz;")

  bform <- bf(y ~ x)
  scode <- stancode(bform, data, brmsfamily("cox", "identity"))
  expect_match2(scode, "target += cox_lpdf(Y | mu, bhaz, cbhaz);")

  bform <- bf(y | bhaz(gr = g) ~ x)
  scode <- stancode(bform, data, brmsfamily("cox"))
  expect_match2(scode, "lprior += dirichlet_lpdf(sbhaz[k] | con_sbhaz[k]);")
  expect_match2(scode, "bhaz[n] = Zbhaz[n] * sbhaz[Jgrbhaz[n]];")
})

test_that("offsets appear in the Stan code", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  scode <- stancode(y ~ x + offset(c), data)
  expect_match2(scode, "+ offsets;")
  scode <- stancode(bf(y ~ a, a ~ offset(log(c + 1)), nl = TRUE),
                         data, prior = prior(normal(0,1), nlpar = a))
  expect_match2(scode, "+ offsets_a;")
})

test_that("prior only models are correctly checked", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  prior <- prior(normal(0, 5), b) + prior("", Intercept)
  expect_error(stancode(y ~ x, data, prior = prior,
                             sample_prior = "only"),
               "Sampling from priors is not possible")
  prior <- prior(normal(0, 5), b) + prior(normal(0, 10), Intercept)
  scode <- stancode(y ~ x, data, prior = prior,
                         sample_prior = "only")
  expect_match2(scode, "lprior += normal_lpdf(Intercept | 0, 10)")
})

test_that("Stan code of mixture model is correct", {
  data <- data.frame(y = 1:10, x = rnorm(10), c = 1)
  data$z <- abs(data$y)

  scode <- stancode(
    bf(y ~ x,  sigma2 ~ x), data,
    family = mixture(gaussian, gaussian),
    sample_prior = TRUE
  )
  expect_match2(scode, "ordered[2] ordered_Intercept;")
  expect_match2(scode, "Intercept_mu2 = ordered_Intercept[2];")
  expect_match2(scode, "lprior += dirichlet_lpdf(theta | con_theta);")
  expect_match2(scode, "ps[1] = log(theta1) + normal_lpdf(Y[n] | mu1[n], sigma1);")
  expect_match2(scode, "ps[2] = log(theta2) + normal_lpdf(Y[n] | mu2[n], sigma2[n]);")
  expect_match2(scode, "target += log_sum_exp(ps);")
  expect_match2(scode, "simplex[2] prior_theta = dirichlet_rng(con_theta);")

  scode <- stancode(bf(z | weights(c) ~ x, shape1 ~ x, theta1 = 1, theta2 = 2),
                         data = data, mixture(Gamma("log"), weibull))
  expect_match(scode, "data \\{[^\\}]*real<lower=0,upper=1> theta1;")
  expect_match(scode, "data \\{[^\\}]*real<lower=0,upper=1> theta2;")
  expect_match2(scode, "ps[1] = log(theta1) + gamma_lpdf(Y[n] | shape1[n], shape1[n] ./ mu1[n]);")
  expect_match2(scode, "target += weights[n] * log_sum_exp(ps);")

  scode <- stancode(bf(abs(y) | se(c) ~ x), data = data,
                         mixture(gaussian, student))
  expect_match2(scode, "ps[1] = log(theta1) + normal_lpdf(Y[n] | mu1[n], se[n]);")
  expect_match2(scode, "ps[2] = log(theta2) + student_t_lpdf(Y[n] | nu2, mu2[n], se[n]);")

  fam <- mixture(gaussian, student, exgaussian)
  scode <- stancode(bf(y ~ x), data = data, family = fam)
  expect_match(scode, "parameters \\{[^\\}]*real Intercept_mu3;")
  expect_match2(scode,
    "ps[2] = log(theta2) + student_t_lpdf(Y[n] | nu2, mu2[n], sigma2);"
  )
  expect_match2(scode,
    "ps[3] = log(theta3) + exp_mod_normal_lpdf(Y[n] | mu3[n] - beta3, sigma3, inv(beta3));"
  )

  scode <- stancode(bf(y ~ x, theta1 ~ x, theta3 ~ x),
                         data = data, family = fam)
  expect_match2(scode, "log_sum_exp_theta = log(exp(theta1[n]) + exp(theta2[n]) + exp(theta3[n]));")
  expect_match2(scode, "theta2 = rep_vector(0.0, N);")
  expect_match2(scode, "theta3[n] = theta3[n] - log_sum_exp_theta;")
  expect_match2(scode, "ps[1] = theta1[n] + normal_lpdf(Y[n] | mu1[n], sigma1);")

  fam <- mixture(cumulative, sratio)
  scode <- stancode(y ~ x, data, family = fam)
  expect_match2(scode, "ordered_logistic_lpmf(Y[n] | mu1[n], Intercept_mu1);")
  expect_match2(scode, "sratio_logit_lpmf(Y[n] | mu2[n], disc2, Intercept_mu2);")

  # censored mixture model
  fam <- mixture(gaussian, gaussian)
  scode <- stancode(y | cens(2, y2 = 2) ~ x, data, fam)
  expect_match2(scode,
    "ps[2] = log(theta2) + normal_lccdf(Y[n] | mu2[n], sigma2);"
  )
  expect_match2(scode, paste0(
    "ps[2] = log(theta2) + log_diff_exp(\n",
    "          normal_lcdf(rcens[n] | mu2[n], sigma2),"
  ))

  # truncated mixture model
  scode <- stancode(y | trunc(3) ~ x, data, fam)
  expect_match2(scode, paste0(
    "ps[1] = log(theta1) + normal_lpdf(Y[n] | mu1[n], sigma1) -\n",
    "        normal_lccdf(lb[n] | mu1[n], sigma1);"
  ))

  # non-linear mixture model
  bform <- bf(y ~ 1) +
    nlf(mu1 ~ eta^2) +
    nlf(mu2 ~ log(eta) + a) +
    lf(eta + a ~ x) +
    mixture(gaussian, nmix = 2)
  bprior <- prior(normal(0, 1), nlpar = "eta") +
    prior(normal(0, 1), nlpar = "a")
  scode <- stancode(bform, data = data, prior = bprior)
  expect_match2(scode, "mu1[n] = (nlp_eta[n] ^ 2);")
  expect_match2(scode, "mu2[n] = (log(nlp_eta[n]) + nlp_a[n]);")
})

test_that("sparse matrix multiplication is applied correctly", {
  data <- data.frame(y = rnorm(10), x = rnorm(10))

  # linear model
  scode <- stancode(
    bf(y ~ x, sparse = TRUE) + lf(sigma ~ x, sparse = TRUE),
    data, prior = prior(normal(0, 5), coef = "Intercept")
  )
  expect_match2(scode, "wX = csr_extract_w(X);")
  expect_match2(scode,
    "mu += csr_matrix_times_vector(rows(X), cols(X), wX, vX, uX, b);"
  )
  expect_match2(scode,
    "uX_sigma[size(csr_extract_u(X_sigma))] = csr_extract_u(X_sigma);"
  )
  expect_match2(scode,
    paste0(
      "sigma += csr_matrix_times_vector(rows(X_sigma), cols(X_sigma), ",
      "wX_sigma, vX_sigma, uX_sigma, b_sigma);"
    )
  )
  expect_match2(scode, "lprior += normal_lpdf(b[1] | 0, 5);")
  expect_match2(scode, "target += normal_lpdf(Y | mu, sigma);")

  # non-linear model
  scode <- stancode(
    bf(y ~ a, lf(a ~ x, sparse = TRUE), nl = TRUE),
    data, prior = prior(normal(0, 1), nlpar = a)
  )
  expect_match2(scode,
    "vX_a[size(csr_extract_v(X_a))] = csr_extract_v(X_a);"
  )
  expect_match2(scode,
    "nlp_a += csr_matrix_times_vector(rows(X_a), cols(X_a), wX_a, vX_a, uX_a, b_a);"
  )
})

test_that("QR decomposition is included in the Stan code", {
  data <- data.frame(y = rnorm(10), x1 = rnorm(10), x2 = rnorm(10))
  bform <- bf(y ~ x1 + x2, decomp = "QR") +
    lf(sigma ~ 0 + x1 + x2, decomp = "QR")

  # simple priors
  scode <- stancode(bform, data, prior = prior(normal(0, 2)))
  expect_match2(scode, "XQ = qr_thin_Q(Xc) * sqrt(N - 1);")
  expect_match2(scode, "b = XR_inv * bQ;")
  expect_match2(scode, "lprior += normal_lpdf(bQ | 0, 2);")
  expect_match2(scode, "XQ * bQ")
  expect_match2(scode, "XR_sigma = qr_thin_R(X_sigma) / sqrt(N - 1);")

  # horseshoe prior
  scode <- stancode(bform, data, prior = prior(horseshoe(1)))
  expect_match2(scode, "target += std_normal_lpdf(zb);")
  expect_match2(scode, "scales = scales_horseshoe(")
  expect_match2(scode, "sdb = scales[(1):(Kc)];")
  expect_match2(scode, "bQ = zb .* sdb;")
})

test_that("Stan code for Gaussian processes is correct", {
  set.seed(1234)
  dat <- data.frame(y = rnorm(40), x1 = rnorm(40), x2 = rnorm(40),
                    z = factor(rep(3:6, each = 10)))

  prior <- prior(gamma(0.1, 0.1), sdgp) +
    prior(gamma(4, 2), sdgp, coef = gpx2x1)
  scode <- stancode(y ~ gp(x1, cov = "matern32") + gp(x2, by = x1, gr = FALSE),
                         dat, prior = prior)
  expect_match2(scode, "lprior += inv_gamma_lpdf(lscale_1[1]")
  expect_match2(scode, "lprior += gamma_lpdf(sdgp_1 | 0.1, 0.1)")
  expect_match2(scode, "lprior += gamma_lpdf(sdgp_2 | 4, 2)")
  expect_match2(scode, "gp_pred_1 = gp_matern32(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);")
  expect_match2(scode, "gp_pred_2 = gp_exp_quad(Xgp_2, sdgp_2[1], lscale_2[1], zgp_2);")
  expect_match2(scode, "Cgp_2 .* gp_pred_2;")

  prior <- prior + prior(normal(0, 1), lscale, coef = gpx1)
  scode <- stancode(y ~ gp(x1, cov = "matern52") + gp(x2, by = x1, cov = "exponential"),
                         data = dat, prior = prior)
  expect_match2(scode, "lprior += normal_lpdf(lscale_1[1][1] | 0, 1)")
  expect_match2(scode, "gp_pred_1 = gp_matern52(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1)")
  expect_match2(scode, "gp_pred_2 = gp_exponential(Xgp_2, sdgp_2[1], lscale_2[1], zgp_2);")
  expect_match2(scode, "+ Cgp_2 .* gp_pred_2[Jgp_2]")

  # non-isotropic GP
  scode <- stancode(y ~ gp(x1, x2, by = z, iso = FALSE), data = dat)
  expect_match2(scode, "lprior += inv_gamma_lpdf(lscale_1[1][2]")
  expect_match2(scode, "lprior += inv_gamma_lpdf(lscale_1[4][2]")

  scode <- stancode(y ~ gp(x1, x2) + gp(x1, by = z, gr = FALSE), data = dat)
  expect_match2(scode, "gp_exp_quad(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1)")
  expect_match2(scode, "mu[Igp_2_2] += Cgp_2_2 .* gp_pred_2_2;")

  # approximate GPS
  scode <- stancode(
    y ~ gp(x1, k = 10, c = 5/4) +
      gp(x2, by = x1, k = 10, c = 5/4, cov = "matern32"),
    data = dat
  )
  expect_match2(scode, "lprior += inv_gamma_lpdf(lscale_1")
  expect_match2(scode,
    "rgp_1 = sqrt(spd_gp_exp_quad(slambda_1, sdgp_1[1], lscale_1[1])) .* zgp_1;"
  )
  expect_match2(scode,
    "rgp_2 = sqrt(spd_gp_matern32(slambda_2, sdgp_2[1], lscale_2[1])) .* zgp_2;"
  )
  expect_match2(scode, "Cgp_2 .* gp_pred_2[Jgp_2]")

  prior <- c(prior(normal(0, 10), lscale, coef = gpx1, nlpar = a),
             prior(gamma(0.1, 0.1), sdgp, nlpar = a),
             prior(normal(0, 1), b, nlpar = a))
  scode <- stancode(bf(y ~ a, a ~ gp(x1), nl = TRUE),
                         data = dat, prior = prior)
  expect_match2(scode, "lprior += normal_lpdf(lscale_a_1[1][1] | 0, 10)")
  expect_match2(scode, "lprior += gamma_lpdf(sdgp_a_1 | 0.1, 0.1)")
  expect_match2(scode, "gp_exp_quad(Xgp_a_1, sdgp_a_1[1], lscale_a_1[1], zgp_a_1)")

  prior <- prior(gamma(2, 2), lscale, coef = gpx1z5, nlpar = "a")
  scode <- stancode(bf(y ~ a, a ~ gp(x1, by = z, gr = TRUE), nl = TRUE),
                         data = dat, prior = prior, silent = TRUE)
  expect_match2(scode,
    "nlp_a[Igp_a_1_1] += Cgp_a_1_1 .* gp_pred_a_1_1[Jgp_a_1_1];"
  )
  expect_match2(scode, "gp_exp_quad(Xgp_a_1_3, sdgp_a_1[3], lscale_a_1[3], zgp_a_1_3)")
  expect_match2(scode, "lprior += gamma_lpdf(lscale_a_1[3][1] | 2, 2);")
  expect_match2(scode, "target += std_normal_lpdf(zgp_a_1_3);")

  # test warnings and errors
  prior <- prior(normal(0, 1), lscale)
  expect_warning(
    stancode(y ~ gp(x1), data = dat, prior = prior),
    "The global prior 'normal(0, 1)' of class 'lscale' will not be used",
    fixed = TRUE
  )
  expect_error(stancode(y ~ gp(x1, cov = "periodic"), data = dat))
})

test_that("Stan code for SAR models is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  W <- matrix(0, nrow = 10, ncol = 10)
  dat2 <- list(W = W)

  scode <- stancode(
    y ~ x + sar(W), data = dat,
    prior = prior(normal(0.5, 1), lagsar),
    data2 = dat2
  )
  expect_match2(scode,
    "target += normal_lagsar_lpdf(Y | mu, sigma, lagsar, Msar, eigenMsar)"
  )
  expect_match2(scode, "lprior += normal_lpdf(lagsar | 0.5, 1)")

  scode <- stancode(
    y ~ x + sar(W, type = "lag"),
    data = dat, family = student(),
    data2 = dat2
  )
  expect_match2(scode,
    "target += student_t_lagsar_lpdf(Y | nu, mu, sigma, lagsar, Msar, eigenMsar)"
  )

  scode <- stancode(y ~ x + sar(W, type = "error"), data = dat,
                         data2 = dat2)
  expect_match2(scode,
    "target += normal_errorsar_lpdf(Y | mu, sigma, errorsar, Msar, eigenMsar)"
  )

  scode <- stancode(
    y ~ x + sar(W, "error"), data = dat, family = student(),
    prior = prior(beta(2, 3), errorsar),
    data2 = dat2
  )
  expect_match2(scode,
    "target += student_t_errorsar_lpdf(Y | nu, mu, sigma, errorsar, Msar, eigenMsar)"
  )
  expect_match2(scode, "lprior += beta_lpdf(errorsar | 2, 3)")

  expect_error(
    stancode(bf(y ~ sar(W), sigma ~ x), data = dat),
    "SAR models are not implemented when predicting 'sigma'"
  )
})

test_that("Stan code for CAR models is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  edges <- cbind(1:10, 10:1)
  W <- matrix(0, nrow = 10, ncol = 10)
  for (i in seq_len(nrow(edges))) {
    W[edges[i, 1], edges[i, 2]] <- 1
  }
  rownames(W) <- seq_len(nrow(W))
  dat2 <- list(W = W)

  scode <- stancode(y ~ x + car(W), dat, data2 = dat2)
  expect_match2(scode, "real<lower=0,upper=1> car;")
  expect_match2(scode, "real sparse_car_lpdf(vector phi")
  expect_match2(scode, "target += sparse_car_lpdf(")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")

  scode <- stancode(y ~ x + car(W, type = "esicar"), dat, data2 = dat2)
  expect_match2(scode, "real sparse_icar_lpdf(vector phi")
  expect_match2(scode, "target += sparse_icar_lpdf(")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")
  expect_match2(scode, "rcar[Nloc] = - sum(zcar)")

  scode <- stancode(y ~ x + car(W, type = "icar"), dat, data2 = dat2)
  expect_match2(scode, "target += -0.5 * dot_self(zcar[edges1] - zcar[edges2])")
  expect_match2(scode, "target += normal_lpdf(sum(zcar) | 0, 0.001 * Nloc)")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")
  expect_match2(scode, "rcar = zcar * sdcar")

  scode <- stancode(y ~ x + car(W, type = "bym2"), dat, data2 = dat2)
  expect_match2(scode, "target += -0.5 * dot_self(zcar[edges1] - zcar[edges2])")
  expect_match2(scode, "target += normal_lpdf(sum(zcar) | 0, 0.001 * Nloc)")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")
  expect_match2(scode, "lprior += beta_lpdf(rhocar | 1, 1)")
  expect_match2(scode, paste0(
    "rcar = (sqrt(1 - rhocar) * nszcar + ",
    "sqrt(rhocar * inv(car_scale)) * zcar) * sdcar"
  ))

  # apply a CAR term on a distributional parameter other than 'mu'
  scode <- stancode(bf(y ~ x, sigma ~ car(W)), dat, data2 = dat2)
  expect_match2(scode, "real sparse_car_lpdf(vector phi")
  expect_match2(scode, "target += sparse_car_lpdf(")
  expect_match2(scode, "sigma[n] += rcar_sigma[Jloc_sigma[n]]")

  # apply shrinkage priors on a CAR term
  scode <- stancode(bf(y ~ x + car(W)), dat, data2 = dat2,
                         prior = prior(horseshoe(main = TRUE), class = b) +
                           prior(horseshoe(), class = sdcar))
  expect_match2(scode, "sdcar = scales[(1+Kc):(Kc+1)][1];")
})

test_that("Stan code for skew_normal models is correct", {
  dat = data.frame(y = rnorm(10), x = rnorm(10))
  scode <- stancode(y ~ x, dat, skew_normal())
  expect_match2(scode, "delta = alpha / sqrt(1 + alpha^2);")
  expect_match2(scode, "omega = sigma / sqrt(1 - sqrt(2 / pi())^2 * delta^2);")
  expect_match2(scode, "mu[n] = mu[n] - omega * delta * sqrt(2 / pi());")

  scode <- stancode(bf(y ~ x, sigma ~ x), dat, skew_normal())
  expect_match2(scode, "omega[n] = sigma[n] / sqrt(1 - sqrt(2 / pi())^2 * delta^2);")
  expect_match2(scode, "mu[n] = mu[n] - omega[n] * delta * sqrt(2 / pi());")

  scode <- stancode(bf(y | se(x) ~ x, alpha ~ x), dat, skew_normal())
  expect_match2(scode, "delta[n] = alpha[n] / sqrt(1 + alpha[n]^2);")
  expect_match2(scode, "omega[n] = se[n] / sqrt(1 - sqrt(2 / pi())^2 * delta[n]^2);")
  expect_match2(scode, "mu[n] = mu[n] - omega[n] * delta[n] * sqrt(2 / pi());")

  scode <- stancode(y ~ x, dat, mixture(skew_normal, nmix = 2))
  expect_match2(scode, "omega1 = sigma1 / sqrt(1 - sqrt(2 / pi())^2 * delta1^2);")
  expect_match2(scode, "mu2[n] = mu2[n] - omega2 * delta2 * sqrt(2 / pi());")
})

test_that("Stan code for missing value terms works correctly", {
  dat = data.frame(y = rnorm(10), x = rnorm(10), g = 1:10, z = 1)
  dat$x[c(1, 3, 9)] <- NA

  bform <- bf(y ~ mi(x)*g) + bf(x | mi() ~ g) + set_rescor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode, "Yl_x[Jmi_x] = Ymi_x;")
  expect_match2(scode, "(bsp_y[1]) * Yl_x[n] + (bsp_y[2]) * Yl_x[n] * Csp_y_1[n];")
  expect_match2(scode, "target += normal_id_glm_lpdf(Yl_x | Xc_x, Intercept_x, b_x, sigma_x);")

  bform <- bf(y ~ mi(x) + (mi(x) | g)) + bf(x | mi() ~ 1) + set_rescor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode,
    "(bsp_y[1] + r_1_y_2[J_1_y[n]]) * Yl_x[n] + r_1_y_1[J_1_y[n]] * Z_1_y_1[n];"
  )

  bform <- bf(y ~ a, a ~ mi(x), nl = TRUE) + bf(x | mi() ~ 1) + set_rescor(FALSE)
  bprior <- prior(normal(0, 1), nlpar = "a", resp = "y")
  scode <- stancode(bform, dat, prior = bprior)
  expect_match2(scode, "nlp_y_a[n] += (bsp_y_a[1]) * Yl_x[n];")
  expect_match2(scode, "lprior += normal_lpdf(bsp_y_a | 0, 1);")

  bform <- bf(y ~ mi(x)*mo(g)) + bf(x | mi() ~ 1) + set_rescor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode, "(bsp_y[3]) * Yl_x[n] * mo(simo_y_2, Xmo_y_2[n]);")

  bform <- bf(y ~ 1, sigma ~ 1) + bf(x | mi() ~ 1) + set_rescor(TRUE)
  scode <- stancode(bform, dat)
  expect_match2(scode, "Yl[n][2] = Yl_x[n];")
  expect_match2(scode, "sigma[n] = transpose([sigma_y[n], sigma_x]);")
  expect_match2(scode, "LSigma[n] = diag_pre_multiply(sigma[n], Lrescor);")

  bform <- bf(x | mi() ~ y, family = "lognormal")
  scode <- stancode(bform, dat)
  expect_match2(scode, "vector<lower=0>[Nmi] Ymi;")

  bform <- bf(y ~ I(log(mi(x))) * g) +
    bf(x | mi() + trunc(lb = 1) ~ y, family = "lognormal")
  scode <- stancode(bform, dat)
  expect_match2(scode, "vector<lower=1>[Nmi_x] Ymi_x;")
  expect_match2(scode,
    "(bsp_y[1]) * (log(Yl_x[n])) + (bsp_y[2]) * (log(Yl_x[n])) * Csp_y_1[n]"
  )

  bform <- bf(y ~ mi(x)*g) +
    bf(x | mi() + cens(z) ~ y, family = "beta")
  scode <- stancode(bform, dat)
  expect_match2(scode, "vector<lower=0,upper=1>[Nmi_x] Ymi_x;")
  expect_match2(scode,
    "target += beta_lpdf(Y_x[Jevent_x[1:Nevent_x]] | mu_x[Jevent_x[1:Nevent_x]] .* phi_x, (1 - mu_x[Jevent_x[1:Nevent_x]]) .* phi_x);"
  )

  # tests #1608
  bform <- bf(y ~ g + mi(x):g + mi(x):mi(z) + mi(z):g) +
    bf(x | mi() ~ 1) +
    bf(z | mi() ~ 1) +
    set_rescor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode,
    "mu_y[n] += (bsp_y[1]) * Yl_x[n] * Csp_y_1[n] + (bsp_y[2]) * Yl_x[n] * Yl_z[n] + (bsp_y[3]) * Yl_z[n] * Csp_y_2[n];"
  )

  bform <- bf(y | mi() ~ mi(x), shape ~ mi(x), family=weibull()) +
    bf(x| mi() ~ z, family=gaussian()) + set_rescor(FALSE)
  scode <- stancode(bform, data = dat)
  expect_match2(scode, "weibull_lpdf(Yl_y | shape_y, mu_y ./ tgamma(1 + 1 ./ shape_y));")
  expect_match2(scode, "shape_y[n] += (bsp_shape_y[1]) * Yl_x[n];")
})

test_that("Stan code for overimputation works correctly", {
  dat = data.frame(y = rnorm(10), x_x = rnorm(10), g = 1:10, z = 1)
  dat$x[c(1, 3, 9)] <- NA
  bform <- bf(y ~ mi(x_x)*g) + bf(x_x | mi(g) ~ 1) + set_rescor(FALSE)
  scode <- stancode(bform, dat, sample_prior = "yes")
  expect_match2(scode, "target += normal_lpdf(Yl_xx | mu_xx, sigma_xx)")
  expect_match2(scode,
    "target += normal_lpdf(Y_xx[Jme_xx] | Yl_xx[Jme_xx], noise_xx[Jme_xx])"
  )
  expect_match2(scode, "vector[N_xx] Yl_xx;")
})

test_that("Missing value terms can be combined with 'subset'", {
  dat <- data.frame(
    y = rnorm(10), x = c(rnorm(9), NA),
    z = rnorm(10), g2 = 10:1,
    g1 = sample(1:5, 10, TRUE),
    s = c(FALSE, rep(TRUE, 9))
  )

  bform <- bf(y ~ mi(x, idx = g1)*mi(z)) +
    bf(x | mi() + index(g2) + subset(s)  ~ 1) +
    bf(z | mi() ~ s) +
    set_rescor(FALSE)
  scode <- stancode(bform, dat)
  expect_match2(scode, "(bsp_y[1]) * Yl_x[idxl_y_x_1[n]]")
  expect_match2(scode, "(bsp_y[2]) * Yl_z[n]")
  expect_match2(scode, "(bsp_y[3]) * Yl_x[idxl_y_x_1[n]] * Yl_z[n]")
  expect_match2(scode, "array[N_y] int idxl_y_x_1;")
})

test_that("Stan code for advanced count data distribution is correct", {
  scode <- stancode(
    count ~ zAge + zBase * Trt + (1|patient),
    data = epilepsy, family = brmsfamily("discrete_weibull")
  )
  expect_match2(scode, "mu = inv_logit(mu);")
  expect_match2(scode, "target += discrete_weibull_lpmf(Y[n] | mu[n], shape);")

  scode <- stancode(
    count ~ zAge + zBase * Trt + (1|patient),
    data = epilepsy, family = brmsfamily("com_poisson")
  )
  expect_match2(scode, "target += com_poisson_log_lpmf(Y[n] | mu[n], shape);")
})

test_that("argument 'stanvars' is handled correctly", {
  bprior <- prior(normal(mean_intercept, 10), class = "Intercept")
  mean_intercept <- 5
  stanvars <- stanvar(mean_intercept)
  scode <- stancode(count ~ Trt, data = epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "real mean_intercept;")

  # define a multi_normal prior with known covariance matrix
  bprior <- prior(multi_normal(M, V), class = "b")
  stanvars <- stanvar(rep(0, 2), "M", scode = "vector[K] M;") +
    stanvar(diag(2), "V", scode = "matrix[K, K] V;")
  scode <- stancode(count ~ Trt + zBase, epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "vector[K] M;")
  expect_match2(scode, "matrix[K, K] V;")

  # define a hierarchical prior on the regression coefficients
  bprior <- set_prior("normal(0, tau)", class = "b") +
    set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)
  stanvars <- stanvar(scode = "real<lower=0> tau;",
                      block = "parameters")
  scode <- stancode(count ~ Trt + zBase, epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "real<lower=0> tau;")
  expect_match2(scode, "lprior += normal_lpdf(b | 0, tau);")

  # ensure that variables are passed to the likelihood of a threaded model
  foo <- 0.5
  stanvars <- stanvar(foo) +
    stanvar(scode = "real<lower=0> tau;",
            block = "parameters", pll_args = "real tau")
  scode <- stancode(count ~ 1, data = epilepsy, family = poisson(),
                         stanvars = stanvars, threads = threading(2),
                         parse = FALSE)
  expect_match2(scode,
    "partial_log_lik_lpmf(array[] int seq, int start, int end, data array[] int Y, real Intercept, data real foo, real tau)"
  )
  expect_match2(scode,
    "reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, Intercept, foo, tau)"
  )

  # specify Stan code in the likelihood part of the model block
  stanvars <- stanvar(scode = "mu += 1.0;", block = "likelihood", position = "start")
  scode <- stancode(count ~ Trt + (1|patient), data = epilepsy,
                         stanvars = stanvars)
  expect_match2(scode, "mu += 1.0;")

  stanvars <- stanvar(scode = "mu += 1.0;", block = "likelihood", position = "start")
  scode <- stancode(count ~ Trt + (1|patient), data = epilepsy,
                         stanvars = stanvars, threads = 2, parse = FALSE)
  expect_match2(scode, "mu += 1.0;")


  # add transformation at the end of a block
  stanvars <- stanvar(scode = "r_1_1 = r_1_1 * 2;",
                      block = "tparameters", position = "end")
  scode <- stancode(count ~ Trt + (1 | patient), epilepsy,
                         stanvars = stanvars)
  expect_match2(scode, "r_1_1 = r_1_1 * 2;\n}")

  # use the non-centered parameterization for 'b'
  # unofficial feature not supported anymore for the time being
  # bprior <- set_prior("target += normal_lpdf(zb | 0, 1)", check = FALSE) +
  #   set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)
  # stanvars <- stanvar(scode = "vector[Kc] zb;", block = "parameters") +
  #   stanvar(scode = "real<lower=0> tau;", block = "parameters") +
  #   stanvar(scode = "vector[Kc] b = zb * tau;",
  #           block="tparameters", name = "b")
  # scode <- stancode(count ~ Trt, epilepsy,
  #                        prior = bprior, stanvars = stanvars)
  # expect_match2(scode, "vector[Kc] b = zb * tau;")

  # stanvars <- stanvar(scode = "vector[Ksp] zbsp;", block = "parameters") +
  #   stanvar(scode = "real<lower=0> tau;", block = "parameters") +
  #   stanvar(scode = "vector[Ksp] bsp = zbsp * tau;",
  #           block = "tparameters", name = "bsp")
  # scode <- stancode(count ~ mo(Base), epilepsy, stanvars = stanvars)
  # expect_match2(scode, "vector[Ksp] bsp = zbsp * tau;")
})

test_that("custom families are handled correctly", {
  dat <- data.frame(size = 10, y = sample(0:10, 20, TRUE), x = rnorm(20))

  # define a custom beta-binomial family
  log_lik_beta_binomial2 <- function(i, prep) {
    mu <- prep$dpars$mu[, i]
    tau <- prep$dpars$tau
    trials <- prep$data$vint1[i]
    y <- prep$data$Y[i]
    beta_binomial2_lpmf(y, mu, tau, trials)
  }
  posterior_predict_beta_binomial2 <- function(i, prep, ...) {
    mu <- prep$dpars$mu[, i]
    tau <- prep$dpars$tau
    trials <- prep$data$vint1[i]
    beta_binomial2_rng(mu, tau, trials)
  }
  posterior_epred_beta_binomial2 <- function(prep) {
    mu <- prep$dpars$mu
    trials <- prep$data$vint1
    trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
    mu * trials
  }
  beta_binomial2 <- custom_family(
    "beta_binomial2",
    dpars = c("mu", "tau"),
    links = c("logit", "log"),
    lb = c(NA, 0),
    type = "int",
    vars = c("vint1[n]", "vreal1[n]"),
    log_lik = log_lik_beta_binomial2,
    posterior_epred = posterior_epred_beta_binomial2,
    posterior_predict = posterior_predict_beta_binomial2
  )
  # define custom stan functions
  # real R is just to also test the vreal addition argument
  stan_funs <- "
    real beta_binomial2_lpmf(int y, real mu, real phi, int N, real R) {
      return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
    }
    int beta_binomial2_rng(real mu, real phi, int N, real R) {
      return beta_binomial_rng(N, mu * phi, (1 - mu) * phi);
    }
  "
  stanvars <- stanvar(scode = stan_funs, block = "functions")
  scode <- stancode(
    y | vint(size) + vreal(size) ~ x, data = dat, family = beta_binomial2,
    prior = prior(gamma(0.1, 0.1), class = "tau"),
    stanvars = stanvars
  )
  expect_match2(scode, "array[N] int vint1;")
  expect_match2(scode, "real<lower=0> tau;")
  expect_match2(scode, "mu = inv_logit(mu);")
  expect_match2(scode, "lprior += gamma_lpdf(tau | 0.1, 0.1);")
  expect_match2(scode,
    "target += beta_binomial2_lpmf(Y[n] | mu[n], tau, vint1[n], vreal1[n]);"
  )

  scode <- stancode(
    bf(y | vint(size) + vreal(size) ~ x, tau ~ x),
    data = dat, family = beta_binomial2, stanvars = stanvars
  )
  expect_match2(scode, "tau = exp(tau);")
  expect_match2(scode,
    "target += beta_binomial2_lpmf(Y[n] | mu[n], tau[n], vint1[n], vreal1[n]);"
  )

  # check custom families in mixture models
  scode <- stancode(
    y | vint(size) + vreal(size) + trials(size) ~ x, data = dat,
    family = mixture(binomial, beta_binomial2),
    stanvars = stanvars
  )
  expect_match2(scode,
    "log(theta2) + beta_binomial2_lpmf(Y[n] | mu2[n], tau2, vint1[n], vreal1[n]);"
  )

  # check custom families in multivariate models
  bform <- bf(
    y | vint(size) + vreal(size) + trials(size) ~ x,
    family = beta_binomial2
  ) + bf(x ~ 1, family = gaussian())
  scode <- stancode(bform, data = dat, stanvars = stanvars)
  expect_match2(scode,
    "target += beta_binomial2_lpmf(Y_y[n] | mu_y[n], tau_y, vint1_y[n], vreal1_y[n]);"
  )

  # check vectorized custom families
  beta_binomial2_vec <- custom_family(
    "beta_binomial2_vec",
    dpars = c("mu", "tau"),
    links = c("logit", "log"),
    lb = c(NA, 0),
    type = "int",
    vars = c("vint1", "vreal1"),
    loop = FALSE
  )
  stan_funs_vec <- "
    real beta_binomial2_vec_lpmf(array[] int y, vector mu, real phi, array[] int N, array[] real R) {
      return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
    }
    int beta_binomial2_rng(real mu, real phi, int N, real R) {
      return beta_binomial_rng(N, mu * phi, (1 - mu) * phi);
    }
  "
  stanvars <- stanvar(scode = stan_funs_vec, block = "functions")
  scode <- stancode(
    y | vint(size) + vreal(size) ~ x, data = dat,
    family = beta_binomial2_vec,
    prior = prior(gamma(0.1, 0.1), class = "tau"),
    stanvars = stanvars
  )
  expect_match2(scode,
    "target += beta_binomial2_vec_lpmf(Y | mu, tau, vint1, vreal1);"
  )
})

test_that("likelihood of distributional beta models is correct", {
  # test issue #404
  dat <- data.frame(prop = rbeta(100, shape1 = 2, shape2 = 2))
  scode <- stancode(
    bf(prop ~ 1, phi ~ 1), data = dat, family = Beta()
  )
  expect_match2(scode, "target += beta_lpdf(Y | mu .* phi, (1 - mu) .* phi);")
})

test_that("student-t group-level effects work without errors", {
  scode <- stancode(count ~ Trt + (1|gr(patient, dist = "st")), epilepsy)
  expect_match2(scode, "dfm_1 = sqrt(df_1 * udf_1);")
  expect_match2(scode, "dfm_1 .* (sd_1[1] * (z_1[1]));")
  expect_match2(scode, "lprior += gamma_lpdf(df_1 | 2, 0.1)")
  expect_match2(scode, "target += inv_chi_square_lpdf(udf_1 | df_1);")

  bprior <- prior(normal(20, 5), class = df, group = patient)
  scode <- stancode(
    count ~ Trt + (Trt|gr(patient, dist = "st")),
    epilepsy, prior = bprior
  )
  expect_match2(scode,
    "r_1 = rep_matrix(dfm_1, M_1) .* scale_r_cor(z_1, sd_1, L_1);"
  )
  expect_match2(scode, "lprior += normal_lpdf(df_1 | 20, 5)")
})

test_that("centering design matrices can be changed correctly", {
  dat <- data.frame(y = 1:10, x = 1:10)
  scode <- stancode(
    bf(y ~ x, center = FALSE), data = dat, family = weibull(),
    prior = prior(normal(0,1), coef = Intercept)
  )
  expect_match2(scode, "mu += X * b;")
  expect_match2(scode, "lprior += normal_lpdf(b[1] | 0, 1);")

  bform <- bf(y ~ eta, nl = TRUE) + lf(eta ~ x, center = TRUE)
  scode <- stancode(bform, data = dat)
  expect_match2(scode, "nlp_eta += Intercept_eta + Xc_eta * b_eta;")
})

test_that("to_vector() is correctly removed from prior of SD parameters", {
  # see https://discourse.mc-stan.org/t/prior-for-sd-generate-parsing-text-error/12292/5
  dat <- data.frame(
    y = rnorm(100),
    ID = 1:10,
    group = rep(1:2, each = 5)
  )
  bform <- bf(
    y ~ 1 + (1 | p | gr(ID, by=group)),
    sigma ~ 1 + (1 | p | gr(ID, by=group))
  )
  bprior <- c(
    prior(normal(0, 0.1), class = sd) ,
    prior(normal(0, 0.01), class = sd, dpar = sigma)
  )
  scode <- stancode(
    bform,
    data = dat,
    prior = bprior,
    sample_prior = TRUE
  )
  expect_match2(scode, "prior_sd_1__1 = normal_rng(0,0.1);")
  expect_match2(scode, "prior_sd_1__2 = normal_rng(0,0.01);")
})

test_that("Dirichlet priors can be flexibly included", {
  # tests issue #1165
  dat <- data.frame(y = rnorm(10), x1 = rnorm(10), x2 = rnorm(10))
  bprior <- prior("dirichlet([1,2]')", class = "b")
  scode <- stancode(y ~ x1 + x2, dat, prior = bprior)
  expect_match2(scode, "simplex[Kc] b;")
})

test_that("threaded Stan code is correct", {
  # tests require cmdstanr which is not yet on CRAN
  skip_on_cran()

  # only run if cmdstan >= 2.29 can be found on the system
  # otherwise the canonicalized code will cause test failures
  # TODO: switch to testing with rstan?
  cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  found_cmdstan <- !brms:::is_try_error(cmdstan_version)
  skip_if_not(found_cmdstan && cmdstan_version >= "2.29.0")
  options(brms.backend = "cmdstanr")

  dat <- data.frame(
    count = rpois(236, lambda = 20),
    visit = rep(1:4, each = 59),
    patient = factor(rep(1:59, 4)),
    Age = rnorm(236),
    Trt = factor(sample(0:1, 236, TRUE)),
    AgeSD = abs(rnorm(236, 1)),
    Exp = sample(1:5, 236, TRUE),
    volume = rnorm(236),
    gender = factor(c(rep("m", 30), rep("f", 29)))
  )

  threads <- threading(2, grainsize = 20)
  bform <- bf(
    count ~ Trt*Age + mo(Exp) + s(Age) + offset(Age) + (1+Trt|visit),
    sigma ~ Trt + gp(Age) + gp(volume, by = Trt)
  )
  scode <- stancode(bform, dat, family = student(), threads = threads)
  expect_match2(scode, "real partial_log_lik_lpmf(array[] int seq, int start,")
  expect_match2(scode, "mu[n] += (bsp[1]) * mo(simo_1, Xmo_1[nn])")
  expect_match2(scode, "ptarget += student_t_lpdf(Y[start:end] | nu, mu, sigma);")
  expect_match2(scode, "+ gp_pred_sigma_1[Jgp_sigma_1[start:end]]")
  expect_match2(scode, ".* gp_pred_sigma_2_1[Jgp_sigma_2_1[which_gp_sigma_2_1]];")
  expect_match2(scode, "sigma[start_at_one(Igp_sigma_2_2[which_gp_sigma_2_2], start)] +=")
  expect_match2(scode, "target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y,")

  scode <- stancode(
    visit ~ cs(Trt) + Age, dat, family = sratio(),
    threads = threads,
  )
  expect_match2(scode, "matrix[N, nthres] mucs = Xcs[start:end] * bcs;")
  expect_match2(scode,
    "ptarget += sratio_logit_lpmf(Y[nn] | mu[n], disc, Intercept")
  expect_match2(scode, " - transpose(mucs[n]));")

  scode <- stancode(
    bf(visit ~ a * Trt ^ b, a ~ mo(Exp), b ~ s(Age), nl = TRUE),
    data = dat, family = Gamma("log"),
    prior = set_prior("normal(0, 1)", nlpar = c("a", "b")),
    threads = threads
  )
  expect_match2(scode, "mu[n] = exp(nlp_a[n] * C_1[nn] ^ nlp_b[n]);")
  expect_match2(scode, "ptarget += gamma_lpdf(Y[start:end] | shape, shape ./ mu);")

  bform <- bf(mvbind(count, Exp) ~ Trt) + set_rescor(TRUE)
  scode <- stancode(bform, dat, gaussian(), threads = threads)
  expect_match2(scode, "ptarget += multi_normal_cholesky_lpdf(Y[start:end] | Mu, LSigma);")

  bform <- bf(brms::mvbind(count, Exp) ~ Trt) + set_rescor(FALSE)
  scode <- stancode(bform, dat, gaussian(), threads = threads)
  expect_match2(scode, "target += reduce_sum(partial_log_lik_count_lpmf, seq_count,")
  expect_match2(scode, "target += reduce_sum(partial_log_lik_Exp_lpmf, seq_Exp,")
  expect_match2(scode,
    "ptarget += normal_id_glm_lpdf(Y_Exp[start:end] | Xc_Exp[start:end], Intercept_Exp, b_Exp, sigma_Exp);"
  )

  scode <- stancode(
    visit ~ Trt, dat, family = mixture(poisson(), nmix = 2),
    threads = threading(4, grainsize = 10, static = TRUE)
  )
  expect_match2(scode, "ps[1] = log(theta1) + poisson_log_lpmf(Y[nn] | mu1[n]);")
  expect_match2(scode, "ptarget += log_sum_exp(ps);")
  expect_match2(scode, "target += reduce_sum_static(partial_log_lik_lpmf,")

  # test that code related to censoring is correct
  scode <- stancode(
    count | cens(Trt) ~ Age, dat, family = lognormal(),
    threads = threading(4)
  )
  expect_match2(scode, "else if (cens[nn] == 1) {")
  expect_match2(scode, "Jrcens[Nrcens] = n;")
  expect_match2(scode,
    "ptarget += lognormal_lcdf(Y[add_int(Jlcens[1:Nlcens], start - 1)] | mu[Jlcens[1:Nlcens]], sigma);"
  )
})

test_that("Un-normalized Stan code is correct", {
  # tests require cmdstanr which is not yet on CRAN
  skip_on_cran()

  # only run if cmdstan >= 2.29 can be found on the system
  # otherwise the canonicalized code will cause test failures
  cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  found_cmdstan <- !brms:::is_try_error(cmdstan_version)
  skip_if_not(found_cmdstan && cmdstan_version >= "2.29.0")
  options(brms.backend = "cmdstanr")

  scode <- stancode(
    count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
    data = epilepsy, family = poisson(),
    prior = prior(student_t(5,0,10), class = b) +
            prior(cauchy(0,2), class = sd),
    normalize = FALSE
  )
  expect_match2(scode, "target += poisson_log_glm_lupmf(Y | Xc, mu, b);")
  expect_match2(scode, "lprior += student_t_lupdf(b | 5, 0, 10);")
  expect_match2(scode, "lprior += student_t_lupdf(Intercept | 3, 1.4, 2.5);")
  expect_match2(scode, "lprior += cauchy_lupdf(sd_1 | 0, 2);")
  expect_match2(scode, "target += std_normal_lupdf(z_1[1]);")

  scode <- stancode(
    count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
    data = epilepsy, family = poisson(),
    prior = prior(student_t(5,0,10), class = b) +
            prior(cauchy(0,2), class = sd),
    normalize = FALSE, threads = threading(2)
  )
  expect_match2(scode, "target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, Xc, b,")
  expect_match2(scode, "Intercept, J_1, Z_1_1, r_1_1, J_2, Z_2_1, r_2_1);")
  expect_match2(scode, "ptarget += poisson_log_glm_lupmf(Y[start:end] | Xc[start:end], mu, b);")
  expect_match2(scode, "lprior += student_t_lupdf(b | 5, 0, 10);")
  expect_match2(scode, "lprior += student_t_lupdf(Intercept | 3, 1.4, 2.5);")
  expect_match2(scode, "lprior += cauchy_lupdf(sd_1 | 0, 2);")
  expect_match2(scode, "target += std_normal_lupdf(z_1[1]);")

  # Check that brms custom distributions stay normalized
  scode <- stancode(
    rating ~ period + carry + cs(treat),
    data = inhaler, family = sratio("cloglog"),
    normalize = FALSE
  )
  expect_match2(scode, "target += sratio_cloglog_lpmf(Y[n] | mu[n], disc, Intercept")
  expect_match2(scode, "- transpose(mucs[n]));")

  # Check that user-specified custom distributions stay normalized
  dat <- data.frame(size = 10, y = sample(0:10, 20, TRUE), x = rnorm(20))

  beta_binomial2 <- custom_family(
      "beta_binomial2",
      dpars = c("mu", "tau"),
      links = c("logit", "log"),
      lb = c(NA, 0),
      type = "int",
      vars = c("vint1[n]", "vreal1[n]"),
  )

  stan_funs <- "
      real beta_binomial2_lpmf(int y, real mu, real phi, int N, real R) {
        return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
      }
    "

  stanvars <- stanvar(scode = stan_funs, block = "functions")

  scode <- stancode(
      y | vint(size) + vreal(size) ~ x, data = dat, family = beta_binomial2,
      prior = prior(gamma(0.1, 0.1), class = "tau"),
      stanvars = stanvars, normalize = FALSE,
  )
  expect_match2(scode, "target += beta_binomial2_lpmf(Y[n] | mu[n], tau, vint1[n], vreal1[n]);")
  expect_match2(scode, "gamma_lupdf(tau | 0.1, 0.1);")
})

# the new array syntax is now used throughout brms
# test_that("Canonicalizing Stan code is correct", {
#   # tests require cmdstanr which is not yet on CRAN
#   skip_on_cran()
#
#   # only run if cmdstan >= 2.29 can be found on the system
#   # otherwise the canonicalized code will cause test failures
#   cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
#   found_cmdstan <- !is_try_error(cmdstan_version)
#   skip_if_not(found_cmdstan && cmdstan_version >= "2.29.0")
#   options(brms.backend = "cmdstanr")
#
#   scode <- stancode(
#     count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
#     data = epilepsy, family = poisson(),
#     prior = prior(student_t(5,0,10), class = b) +
#       prior(cauchy(0,2), class = sd),
#     normalize = FALSE
#   )
#   expect_match2(scode, "array[M_1] vector[N_1] z_1;")
#   expect_match2(scode, "array[M_2] vector[N_2] z_2;")
#
#   model <- "
#   data {
#     int a[5];
#     real b[5];
#     vector[5] c[4];
#   }
#   parameters {
#     real d[5];
#     vector[5] e[4];
#   }
#   "
#   stan_file <- cmdstanr::write_stan_file(model)
#   canonicalized_code <- .canonicalize_stan_model(stan_file, overwrite_file = FALSE)
#   expect_match2(canonicalized_code, "array[5] int a;")
#   expect_match2(canonicalized_code, "array[5] real b;")
#   expect_match2(canonicalized_code, "array[4] vector[5] c;")
#   expect_match2(canonicalized_code, "array[5] real d;")
#   expect_match2(canonicalized_code, "array[4] vector[5] e;")
# })

test_that("Normalizing Stan code works correctly", {
  normalize_stancode <- brms:::normalize_stancode
  expect_equal(
    normalize_stancode("// a\nb;\n  b + c = 4; // kde\ndata"),
    normalize_stancode("// dasflkjldl\n   // adsfadsfa\n b;\n\n  \n  \t\rb + c = 4;\ndata")
  )
  expect_equal(
    normalize_stancode("data /* adfa */ {\nint a;\n /* asdddede \n asdfas \n asf */}\n"),
    normalize_stancode("data {\nint a;\n} /* aa \n adfasdf \n asdfadsf ddd */\n")
  )
  expect_equal(
    normalize_stancode("data \n {\nint a;\n\n }  \t\n"),
    normalize_stancode("data {\nint a;\n} \n")
  )
  expect_equal(
    normalize_stancode("/* \n\n */\na*/"),
    normalize_stancode("a*/")
  )
  expect_equal(
    normalize_stancode("//adsfadf \ra // asdfasdf\r\n"),
    normalize_stancode("a")
  )
  expect_equal(
    normalize_stancode("/* * \n * \n * fg / */hhh"),
    normalize_stancode("hhh")
  )
  expect_equal(
    normalize_stancode("a //b"),
    normalize_stancode("a")
  )
  expect_false(normalize_stancode("// a\ndata {\nint a;\n}\n") ==
                 normalize_stancode("// a\ndata {\nint b;\n}\n"))
  # should not remove single whitespace
  expect_false(normalize_stancode("da ta") ==
                 normalize_stancode("data"))
  # should handle wrong nested comments
  expect_false(normalize_stancode("/* \n\n */\na*/") ==
                 normalize_stancode("b*/"))
})
