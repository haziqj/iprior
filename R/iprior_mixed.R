iprior_mixed <- function(mod, theta0 = NULL, em.maxit = 5, stop.crit = 1e-5,
                         silent = FALSE, control.optim = list()) {
  # Default optim control list -------------------------------------------------
  control.optim_ <- list(
    fnscale = -2,
    trace   = ifelse(isTRUE(silent), 0, 1),
    maxit   = 100,
    REPORT  = 10
  )
  control.optim <- update_control(control.optim, control.optim_)

  # First pass to EM routine ---------------------------------------------------
  cat(paste0("Running ", em.maxit, " initial EM iterations\n"))
  start.time <- Sys.time()
  em.method <- iprior_method_checker(mod, "em")
  if (em.method["em.closed"]) {
    tmp <- iprior_em_closed(mod, em.maxit, stop.crit, silent, theta0,
                            mixed = TRUE)
  }
  if (em.method["em.reg"]) {
    tmp <- iprior_em_reg(mod, em.maxit, stop.crit, silent, theta0)
  }

  # Then pass to direct maximisation routine -----------------------------------
  cat("Now switching to direct optimisation\n")
  res <- iprior_direct(mod, loglik_iprior, tmp$theta, control.optim)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Update time, call, maxit, niter, lb, error, brier --------------------------
  res$time <- time.taken
  res$start.time <- start.time
  res$end.time <- end.time
  res$loglik <- c(tmp$loglik, res$loglik)

  res
}

update_control <- function(arg.list, default.list) {
  default_names <- names(default.list)
  default.list[(arg_names <- names(arg.list))] <- arg.list
  if (length(noNms <- arg_names[!arg_names %in% default_names])) {
    warning("Unknown names in control options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  default.list
}
