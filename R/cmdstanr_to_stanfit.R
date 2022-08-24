# read in stan csvs via cmdstanr and repackage into a stanfit object

read_csv_as_stanfit <- function(files, variables = NULL,
                                sampler_diagnostics = NULL, format = NULL) {

  cmdstanr_read <- cmdstanr::read_cmdstan_csv(files = files,
                                              variables = variables,
                                              sampler_diagnostics = sampler_diagnostics,
                                              format = format)

  n_chains <- nchains(cmdstanr_read$warmup_draws)
  n_iter_warmup <- niterations(cmdstanr_read$warmup_draws)
  n_iter_sample <- niterations(cmdstanr_read$post_warmup_draws)

  # @model_name
  model_name = gsub(".csv", "", basename(files[[1]]))

  # @model_pars
  stanvars <- cmdstanr_read$metadata$stan_variables
  stanvars <- c(stanvars[stanvars != "lp__"], "lp__")

  # @par_dims
  par_dims <- vector("list", length(stanvars))
  names(par_dims) <- stanvars
  par_dims <- lapply(par_dims, function(x) x <- integer(0))

  pdims_num <- sapply(stanvars,
                      function(x) sum(grepl(paste0("^", x, "\\[.*\\]$"), cmdstanr_read$metadata$model_params)))

  par_dims[pdims_num != 0] <- cmdstanr_read$metadata$stan_variable_sizes[stanvars][pdims_num != 0]

  # @mode

  # @sim
  if(n_iter_warmup != 0) {
    cmdstanr_read$warmup_draws <- posterior::as_draws_df(cmdstanr_read$warmup_draws)
    cmdstanr_read$warmup_sampler_diagnostics <-
      posterior::as_draws_df(cmdstanr_read$warmup_sampler_diagnostics)
  }
  if(n_iter_sample != 0) {
    cmdstanr_read$post_warmup_draws <- posterior::as_draws_df(cmdstanr_read$post_warmup_draws)
    cmdstanr_read$post_warmup_sampler_diagnostics <-
      posterior::as_draws_df(cmdstanr_read$post_warmup_sampler_diagnostics)
  }

  draws_full <- rbind(cmdstanr_read$warmup_draws, cmdstanr_read$post_warmup_draws) %>%
    dplyr::relocate("lp__", .after=last_col()) %>%
    as.data.frame

  # manage memory
  cmdstanr_read$warmup_draws <- NULL
  cmdstanr_read$post_warmup_draws <- NULL

  fnames_oi <- colnames(draws_full)[!(colnames(draws_full) %in%
                                        c(".chain", ".iteration", ".draw"))]

  colnames(draws_full) <- gsub("\\[", ".", colnames(draws_full))
  colnames(draws_full) <- gsub("\\]", "", colnames(draws_full))
  colnames(draws_full) <- gsub("\\,", ".", colnames(draws_full))

  samples <- split(draws_full[,!(colnames(draws_full) %in% c(".chain", ".iteration", ".draw"))],
                   draws_full$.chain)
  names(samples) <- NULL

  # manage memory
  rm(draws_full)

  diagnostics_full <- rbind(cmdstanr_read$warmup_sampler_diagnostics,
                            cmdstanr_read$post_warmup_sampler_diagnostics) %>%
    as.data.frame

  # manage memory
  cmdstanr_read$warmup_sampler_diagnostics <- NULL
  cmdstanr_read$post_warmup_sampler_diagnostics <- NULL

  diagnostics_split <- split(diagnostics_full[,!(colnames(diagnostics_full) %in% c(".chain", ".iteration", ".draw"))],
                             diagnostics_full$.chain)
  names(diagnostics_split) <- NULL

  # manage memory
  rm(diagnostics_full)

  # 126-143 lifted from rstan codebase
  values <- c()
  values$algorithm <- cmdstanr_read$metadata$algorithm
  values$engine <- cmdstanr_read$metadata$engine
  values$metric <- cmdstanr_read$metadata$metric

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

  # TODO: add check for existence of post-warmup draws

  #  @sim$sample
  rstan_diagn_order <- c("accept_stat__", "treedepth__", "stepsize__",
                         "divergent__", "n_leapfrog__", "energy__")

  for (i in seq_along(samples)) {
    m <- colMeans(samples[[i]][(n_iter_warmup + 1):(n_iter_warmup + n_iter_sample),])
    rownames(samples[[i]]) <- 1:nrow(samples[[i]])
    attr(samples[[i]], "sampler_params") <- diagnostics_split[[i]][rstan_diagn_order]
    rownames(attr(samples[[i]], "sampler_params")) <- 1:nrow(diagnostics_split[[i]])

    attr(samples[[i]], "adaptation_info") <- paste0("# Step size = ",
                                                    cmdstanr_read$step_size[[i]],
                                                    "\n# Diagonal elements of inverse mass matrix:\n# ",
                                                    paste0(cmdstanr_read$inv_metric[[i]], collapse = ", "),
                                                    "\n# ")

    attr(samples[[i]], "args") <- list(sampler_t = sampler_t,
                                       chain_id = i)

    time_i <- as.double(cmdstanr_read$metadata$time[i,c("warmup", "sampling")])
    names(time_i) <- c("warmup", "sample")
    attr(samples[[i]], "elapsed_time") <- time_i

    attr(samples[[i]], "mean_pars") <- m[-length(m)]
    attr(samples[[i]], "mean_lp__") <- m["lp__"]
  }

  perm_lst <- lapply(1:n_chains, function(id) sample.int(n_iter_sample))

  sim <- list(samples = samples, # some naming to sort out
              iter = cmdstanr_read$metadata$iter_sampling + cmdstanr_read$metadata$iter_warmup,
              thin = cmdstanr_read$metadata$thin,
              warmup = cmdstanr_read$metadata$iter_warmup,
              chains = n_chains,
              n_save = rep(n_iter_sample + n_iter_warmup, n_chains),
              warmup2 = rep(n_iter_warmup, n_chains),
              permutation = perm_lst,
              pars_oi = stanvars,
              dims_oi = par_dims,
              fnames_oi = fnames_oi,
              n_flatnames = length(fnames_oi))

  # @stan_args
  sargs <- list(stan_version_major = as.character(cmdstanr_read$metadata$stan_version_major),
                stan_version_minor = as.character(cmdstanr_read$metadata$stan_version_minor),
                stan_version_patch = as.character(cmdstanr_read$metadata$stan_version_patch),
                model = cmdstanr_read$metadata$model_name,
                start_datetime = gsub(" ", "", cmdstanr_read$metadata$start_datetime),
                method = cmdstanr_read$metadata$method,
                iter = cmdstanr_read$metadata$iter_sampling + cmdstanr_read$metadata$iter_warmup,
                warmup = cmdstanr_read$metadata$iter_warmup,
                save_warmup = cmdstanr_read$metadata$save_warmup,
                thin = cmdstanr_read$metadata$thin,
                engaged = as.character(cmdstanr_read$metadata$adapt_engaged),
                gamma = cmdstanr_read$metadata$gamma,
                delta = cmdstanr_read$metadata$adapt_delta,
                kappa = cmdstanr_read$metadata$kappa,
                t0 = cmdstanr_read$metadata$t0,
                init_buffer = as.character(cmdstanr_read$metadata$init_buffer),
                term_buffer = as.character(cmdstanr_read$metadata$term_buffer),
                window = as.character(cmdstanr_read$metadata$window),
                algorithm = cmdstanr_read$metadata$algorithm,
                engine = cmdstanr_read$metadata$engine,
                max_depth = cmdstanr_read$metadata$max_treedepth,
                metric = cmdstanr_read$metadata$metric,
                metric_file = character(0), # not stored in metadata
                stepsize = NA, # add in loop
                stepsize_jitter = cmdstanr_read$metadata$stepsize_jitter,
                num_chains = as.character(cmdstanr_read$metadata$num_chains),
                chain_id = NA, # add in loop
                file = character(0), # not stored in metadata
                init = NA, # add in loop
                seed = as.character(cmdstanr_read$metadata$seed),
                file = NA, # add in loop
                diagnostic_file = character(0), # not stored in metadata
                refresh = as.character(cmdstanr_read$metadata$refresh),
                sig_figs = as.character(cmdstanr_read$metadata$sig_figs),
                profile_file = cmdstanr_read$metadata$profile_file,
                num_threads = as.character(cmdstanr_read$metadata$threads_per_chain),
                stanc_version = gsub(" ", "", cmdstanr_read$metadata$stanc_version),
                stancflags = character(0), # not stored in metadata
                adaptation_info = NA, # add in loop
                has_time = is.numeric(cmdstanr_read$metadata$time$total),
                time_info = NA, # add in loop
                sampler_t = sampler_t)

  sargs_rep <- replicate(n_chains, sargs, simplify = F)

  for(i in seq_along(sargs_rep)) {
    sargs_rep[[i]]$chain_id <- i
    sargs_rep[[i]]$init <- as.character(cmdstanr_read$metadata$init[i])
    sargs_rep[[i]]$stepsize <- cmdstanr_read$metadata$step_size[i]
    # two 'file' elements: select the second
    file_idx <- which(names(sargs_rep[[i]]) == "file")
    sargs_rep[[i]][[file_idx[2]]] <- files[[i]]

    sargs_rep[[i]]$adaptation_info <- paste0("# Step size = ",
                                             cmdstanr_read$step_size[[i]],
                                             "\n# Diagonal elements of inverse mass matrix:\n# ",
                                             paste0(cmdstanr_read$inv_metric[[i]], collapse = ", "),
                                             "\n# ")

    # note: formatting different from @sim time
    sargs_rep[[i]]$time_info <- paste0(c("#  Elapsed Time: ",
                                         "#                ",
                                         "#                ",
                                         "# "),
                                       c(cmdstanr_read$metadata$time[i, c("warmup", "sampling", "total")], ""),
                                       c(" seconds (Warm-up)", " seconds (Sampling)", " seconds (Total)", ""))
  }

  # @inits: inits are only ever a list
  inits = list()

  # @stanmodel: 262-266 verbatim from read_stan_csv
  null_dso <- new("cxxdso", sig = list(character(0)), dso_saved = FALSE, dso_filename = character(0),
                  modulename = character(0), system = R.version$system, cxxflags = character(0),
                  .CXXDSOMISC = new.env(parent = emptyenv()))
  null_sm <- new("stanmodel", model_name = model_name, model_code = character(0),
                 model_cpp = list(), dso = null_dso)

  # @date
  sdate <- format(file.info(paste0(temp_dir, "/", files[1]))$mtime, "%a %b %d %X %Y")

  new("stanfit",
      model_name = model_name,
      model_pars = stanvars,
      par_dims = par_dims,
      mode = 0L,
      sim = sim,
      inits = inits,
      stan_args = sargs_rep,
      stanmodel = null_sm,
      date = sdate, # not the time of sampling
      .MISC = new.env(parent = emptyenv()))
}
