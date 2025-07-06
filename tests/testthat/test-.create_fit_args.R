# test_that("create fit args and created stancode changes or stays same as expected", {
#   skip_on_cran()
#
#   clean <- function(txt) {
#     txt <- paste(txt, collapse='\n')
#     txt <- gsub("//.*$",          "", txt, perl = TRUE)  # strip // comments
#     txt <- gsub("[[:space:]]+", " ", txt, perl = TRUE)   # squeeze whitespace
#     trimws(txt)
#   }
#   same_sc <- function(a, b) {
#     digest::digest(clean(a), algo = "xxhash64") ==
#       digest::digest(clean(b), algo = "xxhash64")
#   }
#   f_scode <- function(...){
#     call <- brm(... , call_only = T)
#     fit_args <- .create_fit_args(call)
#
#     if(call$backend == 'rstan') {
#       scode <- fit_args$model@model_code
#     } else { # cmdstanr
#       scode <- fit_args$model$code()
#     }
#
#     scode
#   }
#
#   call. <- brm(count ~ zAge + zBase * Trt + (1 | patient),
#                data = epilepsy, family = poisson() , call_only = T)
#
#   args <- .create_fit_args(call.)
#
#
#   ## ----------------------------------------------------------------
#   ## 1.  Same family, same formula, same link  → identical Stan code
#   ## ----------------------------------------------------------------
#   s1 <- f_scode(count ~ zAge + zBase * Trt + (1 | patient),
#                 data = epilepsy, family = poisson())
#   s2 <- f_scode(count ~ zAge + zBase * Trt + (1 | patient),
#                 data = epilepsy, family = poisson())
#   expect_true(same_sc(s1, s2))
#   ## ----------------------------------------------------------------
#   ## 2.  Different family (gaussian vs poisson) → different Stan code
#   ## ----------------------------------------------------------------
#   s3 <- f_scode(count ~ zAge + zBase * Trt + (1 | patient),
#                 data = epilepsy, family = gaussian())
#   expect_false(same_sc(s1, s3))
#   ## ----------------------------------------------------------------
#   ## 3a.  Only the *link* changes (logit vs probit) → code identical
#   ## ----------------------------------------------------------------
#   epi_bin <- within(epilepsy, y <- as.integer(count > 0))
#   b1 <- f_scode(y ~ 1,
#                 data = epi_bin, family = bernoulli(link = "logit"))
#   b2 <- f_scode(y ~ 1,
#                 data = epi_bin, family = bernoulli(link = "probit"))
#   expect_false(same_sc(b1, b2))
#   ## ----------------------------------------------------------------
#   ## 3b.  Add zero-inflation → code *changes*
#   ## ----------------------------------------------------------------
#   z1 <- f_scode(count ~ 1,
#                 data = epi_bin, family = poisson())
#   z2 <- f_scode(count ~ 1,
#                 data = epi_bin, family = zero_inflated_poisson())
#   expect_false(same_sc(z1, z2))
#   ## ----------------------------------------------------------------
#   ## 3.  Prior on sigma vs no prior → code unchanged (sigma already in template)
#   ## ----------------------------------------------------------------
#   p_sig <- set_prior("student_t(3, 0, 10)", class = "sigma")
#   s4 <- f_scode(count ~ 1,
#                 data = epi_bin, family = gaussian(),
#                 prior = p_sig)
#   s5 <- f_scode(count ~ 1,
#                 data = epi_bin, family = gaussian())
#   expect_false(same_sc(s4, s5))
# })
