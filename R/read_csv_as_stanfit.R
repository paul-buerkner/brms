# read in stan CSVs via cmdstanr and repackage into a stanfit object
# efficient replacement of rstan::read_stan_csv
# @param files character vector of CSV files names where draws are stored
# @param variables character vectors of variables to extract draws for
# @param sampler_diagnostics character vectors of diagnostics to extract
# @return an rstan stanfit object
read_csv_as_stanfit <- function(files, variables = NULL,
                                sampler_diagnostics = NULL) {

  csfit <- cmdstanr::read_cmdstan_csv(
    files = files, variables = variables,
    sampler_diagnostics = sampler_diagnostics,
    format = NULL
  )

  # @model_name
  model_name = gsub(".csv", "", basename(files[[1]]))

  # @model_pars
  stanvars <- csfit$metadata$stan_variables
  if ("lp__" %in% stanvars) {
    stanvars <- c(setdiff(stanvars, "lp__"), "lp__")
  }

  if (is.null(variables)) {
    pars_oi <- stanvars
  } else {
    pars_oi <- unique(gsub("\\[.*\\]", "", variables))
  }

  par_names <- csfit$metadata$model_params

  # @par_dims
  par_dims <- vector("list", length(stanvars))

  names(par_dims) <- stanvars
  par_dims <- lapply(par_dims, function(x) x <- integer(0))

  pdims_num <- ulapply(
    stanvars, function(x) sum(grepl(paste0("^", x, "\\[.*\\]$"), par_names))
  )
  par_dims[pdims_num != 0] <-
    csfit$metadata$stan_variable_sizes[stanvars][pdims_num != 0]

  # @mode
  mode <- 0L

  # @sim
  rstan_diagn_order <- c("accept_stat__", "treedepth__", "stepsize__",
                         "divergent__", "n_leapfrog__", "energy__")

  if (!is.null(sampler_diagnostics)) {
    rstan_diagn_order <- rstan_diagn_order[rstan_diagn_order %in% sampler_diagnostics]
  }

  res_vars <- c(".chain", ".iteration", ".draw")
  if ("post_warmup_draws" %in% names(csfit)) {
    # for MCMC samplers
    n_chains <- max(
      nchains(csfit$warmup_draws),
      nchains(csfit$post_warmup_draws)
    )
    n_iter_warmup <- niterations(csfit$warmup_draws)
    n_iter_sample <- niterations(csfit$post_warmup_draws)
    if (n_iter_warmup > 0) {
      csfit$warmup_draws <- as_draws_df(csfit$warmup_draws)
      csfit$warmup_sampler_diagnostics <-
        as_draws_df(csfit$warmup_sampler_diagnostics)
    }
    if (n_iter_sample > 0) {
      csfit$post_warmup_draws <- as_draws_df(csfit$post_warmup_draws)
      csfit$post_warmup_sampler_diagnostics <-
        as_draws_df(csfit$post_warmup_sampler_diagnostics)
    }

    # called 'samples' for consistency with rstan
    samples <- rbind(csfit$warmup_draws, csfit$post_warmup_draws)
    # manage memory
    csfit$warmup_draws <- NULL
    csfit$post_warmup_draws <- NULL

    # prepare sampler diagnostics
    diagnostics <- rbind(csfit$warmup_sampler_diagnostics,
                         csfit$post_warmup_sampler_diagnostics)
    # manage memory
    csfit$warmup_sampler_diagnostics <- NULL
    csfit$post_warmup_sampler_diagnostics <- NULL
    # convert to regular data.frame
    diagnostics <- as.data.frame(diagnostics)
    diag_chain_ids <- diagnostics$.chain
    diagnostics[res_vars] <- NULL

  } else if ("draws" %in% names(csfit)) {
    # for variational inference "samplers"
    n_chains <- 1
    n_iter_warmup <- 0
    n_iter_sample <- niterations(csfit$draws)
    if (n_iter_sample > 0) {
      csfit$draws <- as_draws_df(csfit$draws)
    }

    # called 'samples' for consistency with rstan
    samples <- csfit$draws
    # manage memory
    csfit$draws <- NULL

    # VI has no sampler diagnostics
    diag_chain_ids <- rep(1L, nrow(samples))
    diagnostics <- as.data.frame(matrix(nrow = nrow(samples), ncol = 0))
  }

  # convert to regular data.frame
  samples <- as.data.frame(samples)
  chain_ids <- samples$.chain
  samples[res_vars] <- NULL
  if ("lp__" %in% colnames(samples)) {
    samples <- move2end(samples, "lp__")
  }

  fnames_oi <- colnames(samples)

  colnames(samples) <- gsub("\\[", ".", colnames(samples))
  colnames(samples) <- gsub("\\]", "", colnames(samples))
  colnames(samples) <- gsub("\\,", ".", colnames(samples))

  # split samples into chains
  samples <- split(samples, chain_ids)
  names(samples) <- NULL

  # split diagnostics into chains
  diagnostics <- split(diagnostics, diag_chain_ids)
  names(diagnostics) <- NULL

  #  @sim$sample: largely 113-130 from rstan::read_stan_csv
  values <- list()
  values$algorithm <- csfit$metadata$algorithm
  values$engine <- csfit$metadata$engine
  values$metric <- csfit$metadata$metric

  sampler_t <- NULL
  if (!is.null(values$algorithm)) {
    if (values$algorithm == "rwm" || values$algorithm == "Metropolis") {
      sampler_t <- "Metropolis"
    } else if (values$algorithm == "hmc") {
      if (values$engine == "static") {
        sampler_t <- "HMC"
      } else {
        if (values$metric == "unit_e") {
          sampler_t <- "NUTS(unit_e)"
        } else if (values$metric == "diag_e") {
          sampler_t <- "NUTS(diag_e)"
        } else if (values$metric == "dense_e") {
          sampler_t <- "NUTS(dense_e)"
        }
      }
    }
  }

  adapt_info <- vector("list", 4)
  idx_samples <- (n_iter_warmup + 1):(n_iter_warmup + n_iter_sample)

  for (i in seq_along(samples)) {
    m <- colMeans(samples[[i]][idx_samples, , drop=FALSE])
    rownames(samples[[i]]) <- seq_rows(samples[[i]])
    attr(samples[[i]], "sampler_params") <- diagnostics[[i]][rstan_diagn_order]
    rownames(attr(samples[[i]], "sampler_params")) <- seq_rows(diagnostics[[i]])

    # reformat back to text
    if (is_equal(sampler_t, "NUTS(dense_e)")) {
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

    attr(samples[[i]], "args") <- list(sampler_t = sampler_t, chain_id = i)

    if (NROW(csfit$metadata$time)) {
      time_i <- as.double(csfit$metadata$time[i, c("warmup", "sampling")])
      names(time_i) <- c("warmup", "sample")
      attr(samples[[i]], "elapsed_time") <- time_i
    }

    attr(samples[[i]], "mean_pars") <- m[-length(m)]
    attr(samples[[i]], "mean_lp__") <- m["lp__"]
  }

  perm_lst <- lapply(seq_len(n_chains), function(id) sample.int(n_iter_sample))

  # @sim
  sim <- list(
    samples = samples,
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
    n_flatnames = length(fnames_oi)
  )

  # @stan_args
  sargs <- list(
    stan_version_major = as.character(csfit$metadata$stan_version_major),
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
    sampler_t = sampler_t
  )

  sargs_rep <- replicate(n_chains, sargs, simplify = FALSE)

  for (i in seq_along(sargs_rep)) {
    sargs_rep[[i]]$chain_id <- i
    sargs_rep[[i]]$stepsize <- csfit$metadata$step_size[i]
    sargs_rep[[i]]$init <- as.character(csfit$metadata$init[i])
    # two 'file' elements: select the second
    file_idx <- which(names(sargs_rep[[i]]) == "file")
    sargs_rep[[i]][[file_idx[2]]] <- files[[i]]

    sargs_rep[[i]]$adaptation_info <- adapt_info[[i]]

    if (NROW(csfit$metadata$time)) {
      sargs_rep[[i]]$time_info <- paste0(
        c("#  Elapsed Time: ", "#                ", "#                ", "# "),
        c(csfit$metadata$time[i, c("warmup", "sampling", "total")], ""),
        c(" seconds (Warm-up)", " seconds (Sampling)", " seconds (Total)", "")
      )
    }
  }

  # @stanmodel
  null_dso <- new(
    "cxxdso", sig = list(character(0)), dso_saved = FALSE,
    dso_filename = character(0), modulename = character(0),
    system = R.version$system, cxxflags = character(0),
    .CXXDSOMISC = new.env(parent = emptyenv())
  )
  null_sm <- new(
    "stanmodel", model_name = model_name, model_code = character(0),
    model_cpp = list(), dso = null_dso
  )

  # @date
  sdate <- do.call(max, lapply(files, function(csv) file.info(csv)$mtime))
  sdate <- format(sdate, "%a %b %d %X %Y")

  new(
    "stanfit",
    model_name = model_name,
    model_pars = stanvars,
    par_dims = par_dims,
    mode = mode,
    sim = sim,
    inits = list(),
    stan_args = sargs_rep,
    stanmodel = null_sm,
    date = sdate,  # not the time of sampling
    .MISC = new.env(parent = emptyenv())
  )
}
