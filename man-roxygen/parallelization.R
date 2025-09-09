#' @section Parallelization with multiple CPU cores:
#'
#' \pkg{brms} can make use of multiple CPU cores in parallel to speed
#' up computations in various ways. For efficient use of the
#' available resources it is recommended to only use parallelism to
#' an extend such that the available physical CPUs are not
#' oversubscribed. For example, when you have 8 CPU cores locally
#' available, then you may consider to run 4 chains with 2 threads
#' per chain for best performance if you happen to just run a single
#' model. In case you run a simulation study which requires to run
#' many times a given model, then neither chain nor within-chain
#' parallelization is advisable as the computational resources are
#' already exhausted by the simulation study and any further
#' parallelization beyond the simulation study itself will in fact
#' slow down the overall runtime. Please be aware that for
#' historical reasons the nomenclature of the arguments is possibly
#' confusing. The \code{cores} argument refers to running different
#' chains in parallel and the within-chain parallelization will
#' allocate for each chain as many threads as requested. The
#' requested threads therefore increase the use of overall CPUs in a
#' multiplicative way.
#'
#' For more advanced parallelization (including beyond single model
#' fits), \pkg{brms} also integrates with the \pkg{future}
#' package. Importantly, this enables seamless integration with the
#' \pkg{mirai} parallelization framework through the use of the
#' \pkg{future.mirai} adapter. With \pkg{mirai} local and remote
#' machines can be used in a fully transparent manner to the
#' user. This includes the possibility to use large number of remote
#' machines running in the context of a computer cluster, which are
#' managed with queuing systems. Please refer to the section on
#' distributed computing of
#' \code{\link[mirai:daemons]{mirai::daemons}}.
#'
