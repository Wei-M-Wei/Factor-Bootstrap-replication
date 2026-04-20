# Bootstrap bias correction
bootstrap_bias_correction <- function(Y, X_list, beta_mle, lambda_mle, f_mle,
                                      model_type, num_bootstrap, N, T, R, dimbeta, alpha,  maxiter = 10000, tol = 1e-8, delta, num_cores = 1, use_parallel = FALSE, delta_1 = 0.5) {

  # Fitted values
  z_fit <- lambda_mle %*% t(f_mle)
  for (k in 1:dimbeta) {
    z_fit <- z_fit + X_list[[k]] * beta_mle[k]
  }

  beta_boot <- matrix(NA, dimbeta, num_bootstrap)

  if (num_bootstrap > 0) {
    if (!is.null(num_cores)) {
      n_cores <- min(num_bootstrap, as.integer(num_cores))
    } else {
      n_cores <- min(num_bootstrap, parallel::detectCores(logical = FALSE))
    }
    if (isTRUE(use_parallel) && n_cores > 1) {
      cl <- parallel::makeCluster(n_cores)
      tryCatch({
        parallel::clusterExport(cl, varlist = c("z_fit", "X_list", "dimbeta", "N", "T", "model_type", "R", "maxiter", "tol", "delta"), envir = environment())
        if (exists("NNRPanel_estimate", envir = .GlobalEnv)) {
          parallel::clusterExport(cl, varlist = c("NNRPanel_estimate"), envir = .GlobalEnv)
        }
        parallel::clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1))

        boot_list <- parallel::parLapply(cl, 1:num_bootstrap, function(b) {
          if (model_type == "logit") {
            U <- matrix(runif(N * T), N, T)
            logistic_noise <- log(U / (1 - U))
            Y_boot <- (z_fit + logistic_noise > 0) * 1
          } else if (model_type == "probit") {
            noise <- matrix(rnorm(N * T), N, T)
            Y_boot <- (z_fit + noise > 0) * 1
          }

          data_frame_boot <- matrix(0, N * T, 3 + dimbeta)
          idx <- 1
          for (i in 1:N) {
            for (tt in 1:T) {
              data_frame_boot[idx, 1] <- i
              data_frame_boot[idx, 2] <- tt
              data_frame_boot[idx, 3] <- Y_boot[i, tt]
              for (k in 1:dimbeta) {
                data_frame_boot[idx, 3 + k] <- X_list[[k]][i, tt]
              }
              idx <- idx + 1
            }
          }

          res <- tryCatch({
            result_boot <- NNRPanel_estimate(
              data_frame = data_frame_boot,
              func = model_type,
              delta = delta,
              R_max = R,
              R_true = R,
              R_true_only = TRUE,
              s = 10,
              iter_max = maxiter,
              tol = tol,
              delta_1 = delta_1
            )
            result_boot$beta_fe[1:dimbeta]
          }, error = function(e) {
            rep(NA, dimbeta)
          })
          res
        })

        if (length(boot_list) == num_bootstrap) {
          beta_boot <- do.call(cbind, boot_list)
        }
      }, finally = {
        parallel::stopCluster(cl)
      })
    } else {
      for (b in 1:num_bootstrap) {
        # Generate bootstrap sample
        if (model_type == "logit") {
          U <- matrix(runif(N * T), N, T)
          logistic_noise <- log(U / (1 - U))
          Y_boot <- (z_fit + logistic_noise > 0) * 1
        } else if (model_type == "probit") {
          noise <- matrix(rnorm(N * T), N, T)
          Y_boot <- (z_fit + noise > 0) * 1
        }

        # Create data_frame for bootstrap sample
        data_frame_boot <- matrix(0, N * T, 3 + dimbeta)
        idx <- 1
        for (i in 1:N) {
          for (t in 1:T) {
            data_frame_boot[idx, 1] <- i
            data_frame_boot[idx, 2] <- t
            data_frame_boot[idx, 3] <- Y_boot[i, t]
            for (k in 1:dimbeta) {
              data_frame_boot[idx, 3 + k] <- X_list[[k]][i, t]
            }
            idx <- idx + 1
          }
        }

        # Estimate on bootstrap sample using R_true_only
        tryCatch({
          result_boot <- NNRPanel_estimate(
            data_frame = data_frame_boot,
            func = model_type,
            delta = delta,
            R_max = R,
            R_true = R,
            R_true_only = TRUE,
            s = 10,
            iter_max = maxiter,
            tol = tol
          )
          beta_boot[, b] <- result_boot$beta_fe[1:dimbeta]
        }, error = function(e) {
          beta_boot[, b] <- NA
        })
      }
    }
  }

  # Compute bootstrap bias correction (mean)
  beta_boot_mean <- rowMeans(beta_boot, na.rm = TRUE)
  beta_bc <- 2 * beta_mle - beta_boot_mean

  # Compute bootstrap bias correction (median)
  beta_boot_median <- apply(beta_boot, 1, median, na.rm = TRUE)
  beta_bc_median <- 2 * beta_mle - beta_boot_median

  # Bootstrap CI
  ci_low_1 <- apply(beta_boot - beta_mle, 1, quantile, probs = alpha/2, na.rm = TRUE)
  ci_high_1 <- apply(beta_boot - beta_mle, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE)

  return(list(
    beta_bc = beta_bc,
    beta_bc_median = beta_bc_median,
    beta_boot_raw_mean = beta_boot_mean,
    beta_boot_raw_median = beta_boot_median,
    beta_boot_all = beta_boot,
    ci_low = beta_mle - ci_high_1,
    ci_high = beta_mle - ci_low_1,
    beta_boot = beta_boot
  ))
}

# MLE Estimation
estimate_mle <- function(data, R_true, model_type, maxiter, tol, delta, delta_1 = 0.5) {
  dimbeta <- length(data$X_list)

  result <- tryCatch({
    NNRPanel_estimate(
      data_frame = data$data_frame,
      func = model_type,
      delta = delta,
      R_max = R_true,
      R_true = R_true,
      R_true_only = TRUE,
      s = 100,
      iter_max = maxiter,
      tol = tol,
      delta_1 = delta_1
    )
  }, error = function(e) NULL)

  if (is.null(result)) {
    return(list(
      beta_mle = rep(NA, dimbeta),
      se_beta = rep(NA, dimbeta),
      result = NULL
    ))
  }

  return(list(
    beta_mle = result$beta_fe[1:dimbeta],
    se_beta = if (!is.null(result$std_beta_corr)) result$std_beta_corr[1:dimbeta] else rep(NA, dimbeta),
    result = result
  ))
}

# Split-Panel Jackknife Estimation
estimate_split_pj <- function(mle_result) {
  if (is.null(mle_result$result)) {
    return(list(beta_corr = rep(NA, length(mle_result$beta_mle))))
  }

  dimbeta <- length(mle_result$beta_mle)
  beta_corr <- if (!is.null(mle_result$result$beta_corr_sp)) {
    mle_result$result$beta_corr_sp[1:dimbeta]
  } else {
    rep(NA, dimbeta)
  }

  return(list(beta_corr = beta_corr))
}

# Analytical Bias Correction
estimate_analytical <- function(mle_result) {
  if (is.null(mle_result$result)) {
    return(list(beta_bc_ana = rep(NA, length(mle_result$beta_mle))))
  }

  dimbeta <- length(mle_result$beta_mle)
  beta_bc_ana <- if (!is.null(mle_result$result$beta_corr)) {
    mle_result$result$beta_corr[1:dimbeta]
  } else {
    rep(NA, dimbeta)
  }

  return(list(beta_bc_ana = beta_bc_ana))
}

# Bootstrap Bias Correction
estimate_bootstrap <- function(data, mle_result, R_true, model_type, num_bootstrap, num_cores = NULL, use_parallel = FALSE, alpha, maxiter, tol, delta, delta_1 = 0.5) {
  dimbeta <- length(mle_result$beta_mle)

  # Initialize output
  out <- list(
    beta_bc_boot = rep(NA, dimbeta),
    beta_bc_boot_median = rep(NA, dimbeta),
    beta_boot_raw_mean = rep(NA, dimbeta),
    beta_boot_raw_median = rep(NA, dimbeta),
    beta_boot_all = matrix(NA, nrow = dimbeta, ncol = num_bootstrap),
    ci_low = rep(NA, dimbeta),
    ci_high = rep(NA, dimbeta)
  )

  if (num_bootstrap == 0 || is.null(mle_result$result)) {
    return(out)
  }

  if (is.null(mle_result$result$L_fe) || is.null(mle_result$result$R_fe)) {
    return(out)
  }

  if (!exists("bootstrap_bias_correction")) {
    return(out)
  }

  boot_res <- tryCatch({
    bootstrap_bias_correction(
      Y = data$Y,
      X_list = data$X_list,
      beta_mle = mle_result$beta_mle,
      lambda_mle = mle_result$result$L_fe,
      f_mle = mle_result$result$R_fe,
      model_type = model_type,
      num_bootstrap = num_bootstrap,
      num_cores = num_cores,
      use_parallel = use_parallel,
      N = nrow(data$Y),
      T = ncol(data$Y),
      R = R_true,
      dimbeta = dimbeta,
      alpha = alpha,
      maxiter = maxiter,
      tol = tol,
      delta = delta,
      delta_1 = delta_1
    )
  }, error = function(e) NULL)

  if (!is.null(boot_res)) {
    out$beta_bc_boot <- boot_res$beta_bc
    out$beta_bc_boot_median <- boot_res$beta_bc_median
    out$beta_boot_raw_mean <- boot_res$beta_boot_raw_mean
    out$beta_boot_raw_median <- boot_res$beta_boot_raw_median
    out$beta_boot_all <- boot_res$beta_boot_all
    out$ci_low <- boot_res$ci_low
    out$ci_high <- boot_res$ci_high
  }

  return(out)
}
