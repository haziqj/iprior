iprior_parallel <- function(mod, method = "direct",
                            control = list(silent = FALSE, restarts = TRUE,
                                           no.cores = parallel::detectCores())) {
  if (control$restarts == 1) {
    control$restarts <- parallel::detectCores()
  } else {
    control$no.cores <- min(parallel::detectCores(), control$restarts)
  }
  if (!isTRUE(control$silent)) {
    cat("Performing", control$restarts, "random restarts on", control$no.cores,
        "cores\n")
    snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
    pb <- txtProgressBar(min = 0, max = control$restarts, style = 1)
  } else {
    snow.options.list <- list()
  }

  # The multithreading bit -----------------------------------------------------
  start.time <- Sys.time()
  cl <- parallel::makeCluster(control$no.cores)
  doSNOW::registerDoSNOW(cl)
  res <- foreach::`%dopar%`(
    foreach::foreach(
      i = seq_len(control$restarts),
      .packages = "iprior",
      .options.snow = snow.options.list
    ), {
      new.control          <- control
      new.control$restarts <- 0
      new.control$maxit    <- 3
      new.control$silent   <- TRUE
      tmp <- iprior2(mod, control = new.control, method = method)
      list(
        theta  = tmp$theta,
        loglik = tmp$loglik[length(tmp$loglik)]
      )
    }
  )
  if (!isTRUE(control$silent)) close(pb)
  parallel::stopCluster(cl)

  # Find best starting value ---------------------------------------------------
  tmp <- sapply(res, function(x) x$loglik)
  best.run <- which(tmp == max(tmp))

  # Continue updating the best model -------------------------------------------
  control$restarts <- 0
  control$theta0   <- res[[best.run]]$theta
  control$maxit    <- control$maxit - 3
  res <- iprior2(mod, method = method, control = control)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Update time, call, maxit, niter, lb, error, brier --------------------------
  res$time <- time.taken
  res$start.time <- start.time
  res$end.time <- end.time

  res
}
