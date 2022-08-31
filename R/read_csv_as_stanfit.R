# read in stan csvs via cmdstanr and repackage into a stanfit object

read_csv_as_stanfit <- function(files, variables = NULL, sampler_diagnostics=NULL) {

  csfit <- cmdstanr::read_cmdstan_csv(files = files,
                                      variables = variables,
                                      sampler_diagnostics = sampler_diagnostics,
                                      format = NULL)

  n_chains <- max(nchains(csfit$warmup_draws), nchains(csfit$post_warmup_draws))
  n_iter_warmup <- niterations(csfit$warmup_draws)
  n_iter_sample <- niterations(csfit$post_warmup_draws)

  # @model_name
  model_name = gsub(".csv", "", basename(files[[1]]))

  # @model_pars
  stanvars <- csfit$metadata$stan_variables
  stanvars <- c(stanvars[stanvars != "lp__"], stanvars[stanvars == "lp__"])

  pars_oi <- if(is.null(variables)) {
    stanvars
  } else {
    unique(gsub("\\[.*\\]", "", variables))
  }

  par_names <- csfit$metadata$model_params

  # @par_dims
  par_dims <- vector("list", length(stanvars))

  names(par_dims) <- stanvars
  par_dims <- lapply(par_dims, function(x) x <- integer(0))

  pdims_num <- sapply(stanvars,
                      function(x) sum(grepl(paste0("^", x, "\\[.*\\]$"), par_names)))
  par_dims[pdims_num != 0] <- csfit$metadata$stan_variable_sizes[stanvars][pdims_num != 0]

  # @mode
  mode <- if(n_iter_sample != 0) 0L else 2L

  # @sim
  rstan_diagn_order <- c("accept_stat__", "treedepth__", "stepsize__",
                         "divergent__", "n_leapfrog__", "energy__")

  if(!is.null(sampler_diagnostics)) {
    rstan_diagn_order <- rstan_diagn_order[rstan_diagn_order %in% sampler_diagnostics]
  }

  if(n_iter_warmup != 0) {
    # bind draws
    csfit$warmup_draws <- posterior::as_draws_df(csfit$warmup_draws)
    csfit$warmup_sampler_diagnostics <-
      posterior::as_draws_df(csfit$warmup_sampler_diagnostics)
  }

  if(n_iter_sample != 0) {
    csfit$post_warmup_draws <- posterior::as_draws_df(csfit$post_warmup_draws)
    csfit$post_warmup_sampler_diagnostics <-
      posterior::as_draws_df(csfit$post_warmup_sampler_diagnostics)
  }

  draws_full <- rbind(csfit$warmup_draws, csfit$post_warmup_draws)
  if("lp__" %in% colnames(draws_full)) {
    draws_full <- dplyr::relocate(draws_full, "lp__", .after=last_col())
  }
  draws_full <- as.data.frame(draws_full)

  # manage memory
  csfit$warmup_draws <- NULL
  csfit$post_warmup_draws <- NULL

  fnames_oi <- colnames(draws_full)[!(colnames(draws_full) %in%
                                        c(".chain", ".iteration", ".draw"))]

  colnames(draws_full) <- gsub("\\[", ".", colnames(draws_full))
  colnames(draws_full) <- gsub("\\]", "", colnames(draws_full))
  colnames(draws_full) <- gsub("\\,", ".", colnames(draws_full))

  samples <- split(draws_full[!(colnames(draws_full) %in%
                                   c(".chain", ".iteration", ".draw"))],
                   draws_full$.chain)
  names(samples) <- NULL

  # manage memory
  rm(draws_full)

  # bind diagnostics
  diagnostics_full <- rbind(csfit$warmup_sampler_diagnostics,
                            csfit$post_warmup_sampler_diagnostics)
  diagnostics_full <- as.data.frame(diagnostics_full)

  # manage memory
  csfit$warmup_sampler_diagnostics <- NULL
  csfit$post_warmup_sampler_diagnostics <- NULL

  # split to chains
  diagnostics_split <- split(diagnostics_full[,!(colnames(diagnostics_full) %in%
                                                   c(".chain", ".iteration", ".draw"))],
                             diagnostics_full$.chain)
  names(diagnostics_split) <- NULL

  # manage memory
  rm(diagnostics_full)

  #  @sim$sample: 113-130 verbatim from rstan::read_stan_csv
  values <- c()
  values$algorithm <- csfit$metadata$algorithm
  values$engine <- csfit$metadata$engine
  values$metric <- csfit$metadata$metric

  sampler_t <- NULL
  if (!is.null(values$algorithm) && is.null(values$sampler_t)) {
    if (values$algorithm == "rwm" || values$algorithm ==
        "Metropolis")
      sampler_t <- "Metropolis"
    else if (values$algorithm == "hmc") {
      if (values$engine == "static")
        sampler_t <- "HMC"
      else {
        if (values$metric == "unit_e")
          sampler_t <- "NUTS(unit_e)"
        else if (values$metric == "diag_e")
          sampler_t <- "NUTS(diag_e)"
        else if (values$metric == "dense_e")
          sampler_t <- "NUTS(dense_e)"
      }
    }
  }

  adapt_info <- vector("list", 4)
  idx_samples <- (n_iter_warmup + 1):(n_iter_warmup + n_iter_sample)

  for (i in seq_along(samples)) {
    m <- colMeans(samples[[i]][idx_samples, , drop=FALSE])
    rownames(samples[[i]]) <- seq_len(nrow(samples[[i]]))
    attr(samples[[i]], "sampler_params") <- diagnostics_split[[i]][rstan_diagn_order]
    rownames(attr(samples[[i]], "sampler_params")) <- seq_len(nrow(diagnostics_split[[i]]))

    # reformat back to text
    if(sampler_t == "NUTS(dense_e)") {
      mmatrix_txt <- "\n# Elements of inverse mass matrix:\n# "
      mmat <- paste0(apply(csfit$inv_metric[[i]], 1, paste0, collapse=", "),
                     collapse="\n# ")
    } else {
      mmatrix_txt <- "\n# Diagonal elements of inverse mass matrix:\n# "
      mmat <- paste0(csfit$inv_metric[[i]], collapse = ", ")
    }

    adapt_info[[i]] <- paste0("# Step size = ",
                              csfit$step_size[[i]],
                              mmatrix_txt,
                              mmat, "\n# ")

    attr(samples[[i]], "adaptation_info") <- adapt_info[[i]]

    attr(samples[[i]], "args") <- list(sampler_t = sampler_t,
                                       chain_id = i)

    time_i <- as.double(csfit$metadata$time[i,c("warmup", "sampling")])
    names(time_i) <- c("warmup", "sample")
    attr(samples[[i]], "elapsed_time") <- time_i

    attr(samples[[i]], "mean_pars") <- m[-length(m)]
    attr(samples[[i]], "mean_lp__") <- m["lp__"]
  }

  perm_lst <- lapply(1:n_chains, function(id) sample.int(n_iter_sample))

  # @sim
  sim <- list(samples = samples, # some naming to sort out
              iter = csfit$metadata$iter_sampling + csfit$metadata$iter_warmup,
              thin = csfit$metadata$thin,
              warmup = csfit$metadata$iter_warmup,
              chains = n_chains,
              n_save = rep(n_iter_sample + n_iter_warmup, n_chains),
              warmup2 = rep(n_iter_warmup, n_chains),
              permutation = perm_lst,
              pars_oi = pars_oi,
              dims_oi = par_dims,
              fnames_oi = fnames_oi,
              n_flatnames = length(fnames_oi))

  # @stan_args
  sargs <- list(stan_version_major = as.character(csfit$metadata$stan_version_major),
                stan_version_minor = as.character(csfit$metadata$stan_version_minor),
                stan_version_patch = as.character(csfit$metadata$stan_version_patch),
                model = csfit$metadata$model_name,
                start_datetime = gsub(" ", "", csfit$metadata$start_datetime),
                method = csfit$metadata$method,
                iter = csfit$metadata$iter_sampling + csfit$metadata$iter_warmup,
                warmup = csfit$metadata$iter_warmup,
                save_warmup = csfit$metadata$save_warmup,
                thin = csfit$metadata$thin,
                engaged = as.character(csfit$metadata$adapt_engaged),
                gamma = csfit$metadata$gamma,
                delta = csfit$metadata$adapt_delta,
                kappa = csfit$metadata$kappa,
                t0 = csfit$metadata$t0,
                init_buffer = as.character(csfit$metadata$init_buffer),
                term_buffer = as.character(csfit$metadata$term_buffer),
                window = as.character(csfit$metadata$window),
                algorithm = csfit$metadata$algorithm,
                engine = csfit$metadata$engine,
                max_depth = csfit$metadata$max_treedepth,
                metric = csfit$metadata$metric,
                metric_file = character(0), # not stored in metadata
                stepsize = NA, # add in loop
                stepsize_jitter = csfit$metadata$stepsize_jitter,
                num_chains = as.character(csfit$metadata$num_chains),
                chain_id = NA, # add in loop
                file = character(0), # not stored in metadata
                init = NA, # add in loop
                seed = as.character(csfit$metadata$seed),
                file = NA, # add in loop
                diagnostic_file = character(0), # not stored in metadata
                refresh = as.character(csfit$metadata$refresh),
                sig_figs = as.character(csfit$metadata$sig_figs),
                profile_file = csfit$metadata$profile_file,
                num_threads = as.character(csfit$metadata$threads_per_chain),
                stanc_version = gsub(" ", "", csfit$metadata$stanc_version),
                stancflags = character(0), # not stored in metadata
                adaptation_info = NA, # add in loop
                has_time = is.numeric(csfit$metadata$time$total),
                time_info = NA, # add in loop
                sampler_t = sampler_t)

  sargs_rep <- replicate(n_chains, sargs, simplify = F)

  for(i in seq_along(sargs_rep)) {
    sargs_rep[[i]]$chain_id <- i
    sargs_rep[[i]]$stepsize <- csfit$metadata$step_size[i]
    sargs_rep[[i]]$init <- as.character(csfit$metadata$init[i])
    # two 'file' elements: select the second
    file_idx <- which(names(sargs_rep[[i]]) == "file")
    sargs_rep[[i]][[file_idx[2]]] <- files[[i]]

    sargs_rep[[i]]$adaptation_info <- adapt_info[[i]]

    sargs_rep[[i]]$time_info <- paste0(c("#  Elapsed Time: ",
                                         "#                ",
                                         "#                ",
                                         "# "),
                                       c(csfit$metadata$time[i, c("warmup", "sampling", "total")], ""),
                                       c(" seconds (Warm-up)", " seconds (Sampling)", " seconds (Total)", ""))
  }

  # @stanmodel
  null_dso <- new("cxxdso", sig = list(character(0)), dso_saved = FALSE, dso_filename = character(0),
                  modulename = character(0), system = R.version$system, cxxflags = character(0),
                  .CXXDSOMISC = new.env(parent = emptyenv()))
  null_sm <- new("stanmodel", model_name = model_name, model_code = character(0),
                 model_cpp = list(), dso = null_dso)

  # @date
  sdate <- do.call(max, lapply(files, function(csv) file.info(csv)$mtime))
  sdate <- format(sdate, "%a %b %d %X %Y")

  new("stanfit",
      model_name = model_name,
      model_pars = stanvars,
      par_dims = par_dims,
      mode = mode,
      sim = sim,
      inits = list(),
      stan_args = sargs_rep,
      stanmodel = null_sm,
      date = sdate, # not the time of sampling
      .MISC = new.env(parent = emptyenv()))
}
