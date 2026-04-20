rm(list = ls())
library(NNR)
library(parallel)
library(foreach)
library(doParallel)

# Source bootstrap function if available
source("estimate function.R")


# ============================================================
# PARAMETERS
# ============================================================
set.seed(123)

# Simulation parameters
Tgrid <- c(20,30,40)  # Time periods to test
N <- 30              # Cross-sectional units
S <- 1000            # Number of simulations
num_bootstrap <- 399    # Bootstrap replications (0 = no bootstrap)

# Model parameters
R_true <- 2            # True number of factors
dimbeta <- 1           # Dimension of beta
beta_true <- c( 0.5)    # True beta value
beta_test <- c( 0.7)      # H0 for testing
model_type <- "logit"  # "logit" or "probit"
alpha <- 0.05          # Significance level
maxiter <- 10000       # Maximum iterations for estimation
tol <- 1e-6          # Tolerance for convergence
delta = 1.05

# DGP parameters
r_s <- 1               # Factor loading scale

# DGP Scenario selector
#   1 = Baseline (normal distributions, moderate correlation)
#   2 = High correlation between factors and covariates
#   3 = Outliers in factor loadings
#   4 = Sparse spikes in factors
dgp_scenario <- 3      # Change this to select DGP scenario

# Parallel setup
num_cores <- min(detectCores() - 1, 40)
cat("Setting up parallel cluster with", num_cores, "cores...\n")
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Per-iteration pre-save removed (no longer storing intermediate results)

# ============================================================
# DGP FUNCTION: Generate Data
# ============================================================
# Different scenarios (all with linear factor structure):
#   0: Baseline (normal distributions, moderate correlation)
#   1: Heavy-tailed factors
#   2: High correlation between factors and covariates
#   3: Both heavy-tailed and high correlation
#   4: Outliers in factor loadings
#   5: Heteroskedastic error terms
#   6: Sparse spikes in factors
#   7: Combined challenging scenario
DGP <- function(N, T, R, beta, r_s, model_type, scenario = 1) {
  dimbeta <- length(beta)

  # Set scenario-specific parameters
  if (scenario == 1) {
    # Baseline: normal distributions, moderate correlation
    gamma <- 0.5
    out_frac <- 0
    spike_frac <- 0
    spike_mult <- 1
  } else if (scenario == 2) {
    # High correlation
    gamma <- 0.3
    out_frac <- 0
    spike_frac <- 0
    spike_mult <- 1
  } else if (scenario == 3) {
    # Outliers in factor loadings
    gamma <- 0.5
    out_frac <- 0.1  # 10% of units have outlier loadings
    spike_frac <- 0
    spike_mult <- 1
  } else if (scenario == 4) {
    # Sparse spikes in factors
    gamma <- 0
    out_frac <- 0
    spike_frac <- 0.05  # 5% of time periods have spikes
    spike_mult <- 10    # Spike multiplier
  }

  # Generate factors and loadings
  lambda_i <- matrix(rnorm(N * R), nrow = N, ncol = R)
  f_t <- matrix(rnorm(T * R), nrow = T, ncol = R)

  # Inject outliers in loadings (scenario 4, 7)
  if (out_frac > 0) {
    n_out <- max(1, floor(N * out_frac))
    out_idx <- sample(1:N, size = n_out, replace = FALSE)
    lambda_i[out_idx,] <- matrix(3, nrow = length(out_idx), ncol = R)
  }

  # Inject sparse spikes in factors (scenario 6, 7)
  if (spike_frac > 0) {
    n_spikes <- max(1, floor(T * spike_frac))
    spike_times <- sample(1:T, size = n_spikes, replace = FALSE)
    f_t[spike_times, ] <- f_t[spike_times, ] * spike_mult
  }

  # Interactive fixed effects (linear structure - always linear)
  alpha_true <- r_s * lambda_i %*% t(f_t)

  # Generate covariates (correlated with factors)
  if (scenario == 1 || scenario == 3 || scenario == 4) {
    # Moderate correlation: independent factor structure for X
    lambda_x <- matrix(rnorm(N * R), nrow = N, ncol = R)
    f_x <- matrix(rnorm(T * R), nrow = T, ncol = R)
    Xmat <- lambda_x %*% t(f_x)
  } else {
    # High correlation: X strongly correlated with alpha_true
    Xmat <- gamma * alpha_true
  }

  # Create X array (N x T x dimbeta)
  X_array <- array(0, dim = c(N, T, dimbeta))
  X_list <- list()
  for (k in 1:dimbeta) {
    X_array[, , k] <- Xmat + matrix(rnorm(N * T), nrow = N, ncol = T)
    X_list[[k]] <- X_array[, , k]
  }

  # Compute linear predictor
  Xtheta <- matrix(0, nrow = N, ncol = T)
  for (k in 1:dimbeta) {
    Xtheta <- Xtheta + X_array[, , k] * beta[k]
  }

  # Generate outcome
  if (model_type == "logit") {
    U <- matrix(runif(N * T), nrow = N, ncol = T)
    logistic_noise <- log(U / (1 - U))
    Y <- (Xtheta + alpha_true + logistic_noise > 0) * 1
  } else if (model_type == "probit") {
    noise <- matrix(rnorm(N * T), nrow = N, ncol = T)
    Y <- (Xtheta + alpha_true + noise > 0) * 1
  }

  # Convert to data_frame format [ID, Time, Y, X1, X2, ...]
  data_frame <- matrix(0, N * T, 3 + dimbeta)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:T) {
      data_frame[idx, 1] <- i
      data_frame[idx, 2] <- t
      data_frame[idx, 3] <- Y[i, t]
      for (k in 1:dimbeta) {
        data_frame[idx, 3 + k] <- X_array[i, t, k]
      }
      idx <- idx + 1
    }
  }

  return(list(
    Y = Y,
    X = X_array,
    X_list = X_list,
    data_frame = data_frame,
    lambda_true = lambda_i,
    f_true = f_t,
    alpha_true = alpha_true
  ))
}


# ============================================================
# STORAGE ARRAYS
# ============================================================
num_T <- length(Tgrid)

results <- list(
  beta_mle = array(NA, dim = c(num_T, S, dimbeta)),
  beta_corr = array(NA, dim = c(num_T, S, dimbeta)),
  beta_bc_ana = array(NA, dim = c(num_T, S, dimbeta)),
  beta_bc_boot = array(NA, dim = c(num_T, S, dimbeta)),
  beta_bc_boot_median = array(NA, dim = c(num_T, S, dimbeta)),
  beta_boot_all = array(NA, dim = c(num_T, S, dimbeta, num_bootstrap)),
  ci_low = array(NA, dim = c(num_T, S, dimbeta)),
  ci_high = array(NA, dim = c(num_T, S, dimbeta)),
  Var_beta = array(NA, dim = c(num_T, S, dimbeta)),
  num_factors = matrix(NA, nrow = num_T, ncol = S),
  coverage_mle = matrix(NA, num_T, dimbeta),
  coverage_corr = matrix(NA, num_T, dimbeta),
  coverage_ana = matrix(NA, num_T, dimbeta),
  coverage_boot = matrix(NA, num_T, dimbeta),
  coverage_boot_median = matrix(NA, num_T, dimbeta),
  reject_mle = matrix(NA, num_T, dimbeta),
  reject_corr = matrix(NA, num_T, dimbeta),
  reject_ana = matrix(NA, num_T, dimbeta),
  reject_boot = matrix(NA, num_T, dimbeta),
  reject_boot_median = matrix(NA, num_T, dimbeta)
)

# ============================================================
# SIMULATION LOOP
# ============================================================
cat("\n========================================\n")
cat("MONTE CARLO SIMULATION\n")
cat("========================================\n")
cat("N:", N, "\n")
cat("T grid:", Tgrid, "\n")
cat("Simulations (S):", S, "\n")
cat("True R:", R_true, "\n")
cat("True beta:", beta_true, "\n")
cat("Model:", model_type, "\n")
cat("Bootstrap replications:", num_bootstrap, "\n")
cat("DGP Scenario:", dgp_scenario, "\n")


for (model_type  in c('logit', 'probit')) {
  for(dgp_scenario in c(1,2,3,4)){
    for (t_idx in 1:num_T){
      T_val <- Tgrid[t_idx]
      cat(sprintf("================ Processing T = %d ================\n", T_val))

      start_time <- Sys.time()

      # Run simulations in parallel (direct foreach assignment)
      sim_results <- foreach(
        jsim = 1:S,
        .packages = c("NNR", "MASS"),
        .errorhandling = 'pass',
        .export = c("DGP", "estimate_mle", "estimate_split_pj", "estimate_analytical",
                    "estimate_bootstrap",
                    "R_true", "beta_true", "r_s", "model_type",
                    "num_bootstrap", "alpha", "maxiter", "tol", "dgp_scenario")
      ) %dopar% {
        set.seed(jsim)
        # Generate data with selected DGP scenario
        data <- DGP(N = N, T = T_val, R = R_true, beta = beta_true,
                    r_s = r_s, model_type = model_type, scenario = dgp_scenario)

        # Apply different estimation methods separately
        # 1. MLE
        mle_result <- estimate_mle(data = data, R_true = R_true,
                                   model_type = model_type,
                                   maxiter = maxiter, tol = tol, delta = delta)
        # Safely capture estimated number of factors (if returned by estimator)
        nf_est <- NA
        if (!is.null(mle_result$result) && !is.null(mle_result$result$num_factor_est)) {
          nf_est <- mle_result$result$num_factor_est
        }
        # 2. Split-panel Jackknife
        split_pj_result <- estimate_split_pj(mle_result = mle_result)

        # 3. Analytical Bias Correction
        analytical_result <- estimate_analytical(mle_result = mle_result)

        # 4. Bootstrap Bias Correction
        bootstrap_result <- estimate_bootstrap(data = data, mle_result = mle_result,
                                               R_true = R_true, model_type = model_type,
                                               num_bootstrap = num_bootstrap, alpha = alpha, maxiter = maxiter, tol = tol, delta = delta)

        # Combine all results
        out <- list(
          beta_mle = mle_result$beta_mle,
          beta_corr = split_pj_result$beta_corr,
          beta_bc_ana = analytical_result$beta_bc_ana,
          beta_bc_boot = bootstrap_result$beta_bc_boot,
          beta_bc_boot_median = bootstrap_result$beta_bc_boot_median,
          beta_boot_raw_mean = bootstrap_result$beta_boot_raw_mean,
          beta_boot_raw_median = bootstrap_result$beta_boot_raw_median,
          beta_boot_all = bootstrap_result$beta_boot_all,
          se_beta = mle_result$se_beta,
          ci_low = bootstrap_result$ci_low,
          ci_high = bootstrap_result$ci_high,
          num_factor_est = nf_est
        )

        out
      }

      end_time <- Sys.time()
      elapsed <- difftime(end_time, start_time, units = "secs")
      cat(sprintf("  Completed %d simulations in %.1f seconds\n", length(sim_results), as.numeric(elapsed)))

      # Extract and store results with safety checks
      safe_extract <- function(x, expected_len) {
        if (is.null(x)) return(rep(NA_real_, expected_len))
        if (is.numeric(x)) {
          x_len <- length(x)
          if (x_len == expected_len) {
            return(as.numeric(x))
          } else if (x_len > expected_len) {
            return(as.numeric(x)[1:expected_len])
          } else {
            # Pad with NA if too short
            result <- rep(NA_real_, expected_len)
            result[1:x_len] <- as.numeric(x)
            return(result)
          }
        }
        return(rep(NA_real_, expected_len))
      }

      error_count <- 0
      for (jsim in 1:S) {
        # Skip error results or malformed entries
        if (inherits(sim_results[[jsim]], "error") || is.null(sim_results[[jsim]]) || !is.list(sim_results[[jsim]])) {
          error_count <- error_count + 1
          next
        }
        # Safely extract each component with proper length checking
        results$beta_mle[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$beta_mle, dimbeta)
        results$beta_corr[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$beta_corr, dimbeta)
        results$beta_bc_ana[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$beta_bc_ana, dimbeta)
        results$beta_bc_boot[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$beta_bc_boot, dimbeta)
        results$beta_bc_boot_median[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$beta_bc_boot_median, dimbeta)
        results$num_factors[t_idx, jsim] <- safe_extract(sim_results[[jsim]]$num_factor_est, 1)
        if (num_bootstrap > 0) {
          boot_all <- sim_results[[jsim]]$beta_boot_all
          if (!is.null(boot_all) && is.matrix(boot_all) && nrow(boot_all) == dimbeta && ncol(boot_all) == num_bootstrap) {
            results$beta_boot_all[t_idx, jsim, , ] <- boot_all
          }
        }
        se_beta <- safe_extract(sim_results[[jsim]]$se_beta, dimbeta)
        results$Var_beta[t_idx, jsim, ] <- se_beta^2
        # Save returned bootstrap CIs (if provided by estimate_bootstrap)
        results$ci_low[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$ci_low, dimbeta)
        results$ci_high[t_idx, jsim, ] <- safe_extract(sim_results[[jsim]]$ci_high, dimbeta)
      }
      if (error_count > 0) cat(sprintf("  Warning: %d simulations had errors or invalid results\n", error_count))

      # Compute coverage and rejection rates
      z_alpha <- qnorm(1 - alpha/2)

      for (j in 1:dimbeta) {
        se_j <- sqrt(results$Var_beta[t_idx, , j])

        # MLE
        ci_low <- results$beta_mle[t_idx, , j] - z_alpha * se_j
        ci_high <- results$beta_mle[t_idx, , j] + z_alpha * se_j
        results$coverage_mle[t_idx, j] <- mean(beta_true[j] >= ci_low & beta_true[j] <= ci_high, na.rm = TRUE)
        results$reject_mle[t_idx, j] <- mean(beta_test[j] < ci_low | beta_test[j] > ci_high, na.rm = TRUE)

        # Split-panel Jackknife
        ci_low <- results$beta_corr[t_idx, , j] - z_alpha * se_j
        ci_high <- results$beta_corr[t_idx, , j] + z_alpha * se_j
        results$coverage_corr[t_idx, j] <- mean(beta_true[j] >= ci_low & beta_true[j] <= ci_high, na.rm = TRUE)
        results$reject_corr[t_idx, j] <- mean(beta_test[j] < ci_low | beta_test[j] > ci_high, na.rm = TRUE)

        # Analytical BC
        ci_low <- results$beta_bc_ana[t_idx, , j] - z_alpha * se_j
        ci_high <- results$beta_bc_ana[t_idx, , j] + z_alpha * se_j
        results$coverage_ana[t_idx, j] <- mean(beta_true[j] >= ci_low & beta_true[j] <= ci_high, na.rm = TRUE)
        results$reject_ana[t_idx, j] <- mean(beta_test[j] < ci_low | beta_test[j] > ci_high, na.rm = TRUE)

        # Bootstrap BC: use CIs returned by the bootstrap function when available
        if (num_bootstrap > 0) {
          ci_low_ret <- results$ci_low[t_idx, , j]
          ci_high_ret <- results$ci_high[t_idx, , j]
          results$coverage_boot[t_idx, j] <- mean(beta_true[j] >= ci_low_ret & beta_true[j] <= ci_high_ret, na.rm = TRUE)
          results$reject_boot[t_idx, j] <- mean(beta_test[j] < ci_low_ret | beta_test[j] > ci_high_ret, na.rm = TRUE)

          # For median-based bootstrap BC use the returned CIs as well (if bootstrap function
          # provided separate median CIs you can store/use them; otherwise we reuse returned CIs).
          ci_low_med <- results$ci_low[t_idx, , j]
          ci_high_med <- results$ci_high[t_idx, , j]
          results$coverage_boot_median[t_idx, j] <- mean(beta_true[j] >= ci_low_med & beta_true[j] <= ci_high_med, na.rm = TRUE)
          results$reject_boot_median[t_idx, j] <- mean(beta_test[j] < ci_low_med | beta_test[j] > ci_high_med, na.rm = TRUE)
        }
      }

      # Print results for this T immediately
      cat("\n")
      cat("================================================================================\n")
      cat(sprintf("Results for T = %d\n", T_val))
      cat("================================================================================\n")

      for (j in 1:dimbeta) {
        cat(sprintf("\n--- Covariate %d (True beta = %.3f) ---\n", j, beta_true[j]))
        cat("Method     \tBias\t\tStd\t\tCR\t\tRej\n")
        cat("--------------------------------------------------------------------------------\n")

        # MLE
        beta_mle <- results$beta_mle[t_idx, , j]
        mean_mle <- mean(beta_mle, na.rm = TRUE)
        std_mle <- sd(beta_mle, na.rm = TRUE)
        bias_mle <- (mean_mle - beta_true[j]) / abs(beta_true[j])
        cat(sprintf("MLE       \t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                    bias_mle, std_mle, results$coverage_mle[t_idx, j], results$reject_mle[t_idx, j]))

        # Split-panel Jackknife
        beta_corr <- results$beta_corr[t_idx, , j]
        mean_corr <- mean(beta_corr, na.rm = TRUE)
        std_corr <- sd(beta_corr, na.rm = TRUE)
        bias_corr <- (mean_corr - beta_true[j]) / abs(beta_true[j])
        cat(sprintf("SplitPJ   \t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                    bias_corr, std_corr, results$coverage_corr[t_idx, j], results$reject_corr[t_idx, j]))

        # Analytical BC
        beta_ana <- results$beta_bc_ana[t_idx, , j]
        mean_ana <- mean(beta_ana, na.rm = TRUE)
        std_ana <- sd(beta_ana, na.rm = TRUE)
        bias_ana <- (mean_ana - beta_true[j]) / abs(beta_true[j])
        cat(sprintf("Analytical\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                    bias_ana, std_ana, results$coverage_ana[t_idx, j], results$reject_ana[t_idx, j]))

        # Bootstrap BC (mean)
        if (num_bootstrap > 0) {
          beta_boot <- results$beta_bc_boot[t_idx, , j]
          mean_boot <- mean(beta_boot, na.rm = TRUE)
          std_boot <- sd(beta_boot, na.rm = TRUE)
          bias_boot <- (mean_boot - beta_true[j]) / abs(beta_true[j])
          cat(sprintf("Boot-Mean\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                      bias_boot, std_boot, results$coverage_boot[t_idx, j], results$reject_boot[t_idx, j]))

          # Bootstrap BC (median)
          beta_boot_med <- results$beta_bc_boot_median[t_idx, , j]
          mean_boot_med <- mean(beta_boot_med, na.rm = TRUE)
          std_boot_med <- sd(beta_boot_med, na.rm = TRUE)
          bias_boot_med <- (mean_boot_med - beta_true[j]) / abs(beta_true[j])
          cat(sprintf("Boot-Med\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                      bias_boot_med, std_boot_med, results$coverage_boot_median[t_idx, j], results$reject_boot_median[t_idx, j]))
        }
      }
      # Print estimated factor-count summary for this T
      nf_vec <- results$num_factors[t_idx, ]
      mean_nf <- mean(nf_vec, na.rm = TRUE)
      sd_nf <- sd(nf_vec, na.rm = TRUE)
      prop_eq <- if (all(is.na(nf_vec))) NA else mean(nf_vec == R_true, na.rm = TRUE)
      cat(sprintf("Estimated factors (mean ± sd): %.3f ± %.3f; Pr(==R_true=%d)=%.3f\n",
                  mean_nf, sd_nf, R_true, ifelse(is.na(prop_eq), 0, prop_eq)))

      cat("================================================================================\n\n")
    }

    # ============================================================
    # FINAL SUMMARY (All T values)
    # ============================================================
    cat("\n========================================\n")
    cat("FINAL SUMMARY - ALL T VALUES\n")
    cat("========================================\n")
    cat(sprintf("Model: %s, N = %d, S = %d\n", model_type, N, S))
    cat(sprintf("True beta: %s\n", paste(beta_true, collapse = ", ")))
    cat(sprintf("H0: beta = %s\n\n", paste(beta_test, collapse = ", ")))

    cat("================================================================================\n")

    for (j in 1:dimbeta) {
      cat(sprintf("\n--- Covariate %d (True beta = %.3f) ---\n", j, beta_true[j]))
      cat("T\tMethod     \tBias\t\tStd\t\tCR\t\tRej\n")
      cat("--------------------------------------------------------------------------------\n")

      for (t_idx in 1:num_T) {
        T_val <- Tgrid[t_idx]

        # MLE
        beta_mle <- results$beta_mle[t_idx, , j]
        mean_mle <- mean(beta_mle, na.rm = TRUE)
        std_mle <- sd(beta_mle, na.rm = TRUE)
        bias_mle <- (mean_mle - beta_true[j]) / abs(beta_true[j])

        cat(sprintf("%d\tMLE       \t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                    T_val, bias_mle, std_mle, results$coverage_mle[t_idx, j], results$reject_mle[t_idx, j]))

        # Split-panel Jackknife
        beta_corr <- results$beta_corr[t_idx, , j]
        mean_corr <- mean(beta_corr, na.rm = TRUE)
        std_corr <- sd(beta_corr, na.rm = TRUE)
        bias_corr <- (mean_corr - beta_true[j]) / abs(beta_true[j])

        cat(sprintf("\tSplitPJ   \t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                    bias_corr, std_corr, results$coverage_corr[t_idx, j], results$reject_corr[t_idx, j]))

        # Analytical BC
        beta_ana <- results$beta_bc_ana[t_idx, , j]
        mean_ana <- mean(beta_ana, na.rm = TRUE)
        std_ana <- sd(beta_ana, na.rm = TRUE)
        bias_ana <- (mean_ana - beta_true[j]) / abs(beta_true[j])

        cat(sprintf("\tAnalytical\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                    bias_ana, std_ana, results$coverage_ana[t_idx, j], results$reject_ana[t_idx, j]))

        # Bootstrap BC (mean)
        if (num_bootstrap > 0) {
          beta_boot <- results$beta_bc_boot[t_idx, , j]
          mean_boot <- mean(beta_boot, na.rm = TRUE)
          std_boot <- sd(beta_boot, na.rm = TRUE)
          bias_boot <- (mean_boot - beta_true[j]) / abs(beta_true[j])

          cat(sprintf("\tBoot-Mean\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                      bias_boot, std_boot, results$coverage_boot[t_idx, j], results$reject_boot[t_idx, j]))

          # Bootstrap BC (median)
          beta_boot_med <- results$beta_bc_boot_median[t_idx, , j]
          mean_boot_med <- mean(beta_boot_med, na.rm = TRUE)
          std_boot_med <- sd(beta_boot_med, na.rm = TRUE)
          bias_boot_med <- (mean_boot_med - beta_true[j]) / abs(beta_true[j])

          cat(sprintf("\tBoot-Med\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\n",
                      bias_boot_med, std_boot_med, results$coverage_boot_median[t_idx, j], results$reject_boot_median[t_idx, j]))
        }

        cat("--------------------------------------------------------------------------------\n")
        # Print estimated factor-count summary for this T in the final summary
        nf_vec <- results$num_factors[t_idx, ]
        mean_nf <- mean(nf_vec, na.rm = TRUE)
        sd_nf <- sd(nf_vec, na.rm = TRUE)
        prop_eq <- if (all(is.na(nf_vec))) NA else mean(nf_vec == R_true, na.rm = TRUE)
        cat(sprintf("  Estimated factors (mean ± sd): %.3f ± %.3f; Pr(==R_true=%d)=%.3f\n",
                    mean_nf, sd_nf, R_true, ifelse(is.na(prop_eq), 0, prop_eq)))
      }
    }

    cat("\n*** SIMULATION COMPLETE ***\n")

    # ============================================================
    # SAVE RESULTS
    # ============================================================
    cat("\n========================================\n")
    cat("SAVING RESULTS\n")
    cat("========================================\n")

    # Create summary table
    summary_table <- data.frame(
      T = integer(),
      Covariate = integer(),
      Method = character(),
      Bias = numeric(),
      Std = numeric(),
      Coverage = numeric(),
      Rejection = numeric(),
      stringsAsFactors = FALSE
    )

    for (t_idx in 1:num_T) {
      T_val <- Tgrid[t_idx]

      for (j in 1:dimbeta) {
        # MLE
        beta_mle <- results$beta_mle[t_idx, , j]
        summary_table <- rbind(summary_table, data.frame(
          T = T_val,
          Covariate = j,
          Method = "MLE",
          Bias = (mean(beta_mle, na.rm = TRUE) - beta_true[j]) / abs(beta_true[j]),
          Std = sd(beta_mle, na.rm = TRUE),
          Coverage = results$coverage_mle[t_idx, j],
          Rejection = results$reject_mle[t_idx, j]
        ))

        # Split-panel Jackknife
        beta_corr <- results$beta_corr[t_idx, , j]
        summary_table <- rbind(summary_table, data.frame(
          T = T_val,
          Covariate = j,
          Method = "SplitPJ",
          Bias = (mean(beta_corr, na.rm = TRUE) - beta_true[j]) / abs(beta_true[j]),
          Std = sd(beta_corr, na.rm = TRUE),
          Coverage = results$coverage_corr[t_idx, j],
          Rejection = results$reject_corr[t_idx, j]
        ))

        # Analytical BC
        beta_ana <- results$beta_bc_ana[t_idx, , j]
        summary_table <- rbind(summary_table, data.frame(
          T = T_val,
          Covariate = j,
          Method = "Analytical",
          Bias = (mean(beta_ana, na.rm = TRUE) - beta_true[j]) / abs(beta_true[j]),
          Std = sd(beta_ana, na.rm = TRUE),
          Coverage = results$coverage_ana[t_idx, j],
          Rejection = results$reject_ana[t_idx, j]
        ))

        if (num_bootstrap > 0) {
          # Bootstrap Mean
          beta_boot <- results$beta_bc_boot[t_idx, , j]
          summary_table <- rbind(summary_table, data.frame(
            T = T_val,
            Covariate = j,
            Method = "Boot-Mean",
            Bias = (mean(beta_boot, na.rm = TRUE) - beta_true[j]) / abs(beta_true[j]),
            Std = sd(beta_boot, na.rm = TRUE),
            Coverage = results$coverage_boot[t_idx, j],
            Rejection = results$reject_boot[t_idx, j]
          ))

          # Bootstrap Median
          beta_boot_med <- results$beta_bc_boot_median[t_idx, , j]
          summary_table <- rbind(summary_table, data.frame(
            T = T_val,
            Covariate = j,
            Method = "Boot-Median",
            Bias = (mean(beta_boot_med, na.rm = TRUE) - beta_true[j]) / abs(beta_true[j]),
            Std = sd(beta_boot_med, na.rm = TRUE),
            Coverage = results$coverage_boot_median[t_idx, j],
            Rejection = results$reject_boot_median[t_idx, j]
          ))
        }
      }
    }

    # Save summary table
    filename_summary <- sprintf("%s_summary_N%d_S%d_scenario%d.csv", model_type, N, S, dgp_scenario)
    write.csv(summary_table, filename_summary, row.names = FALSE)
    cat(sprintf("Summary table saved to: %s\n", filename_summary))

    # Save all results as RData
    filename_rdata <- sprintf("%s_results_N%d_S%d.RData", model_type, N, S)
    save(results, file = filename_rdata)
    cat(sprintf("Full results saved to: %s\n", filename_rdata))

    # Save estimated factor count summary
    numfactor_table <- data.frame(
      T = integer(),
      MeanNumFactors = numeric(),
      SdNumFactors = numeric(),
      PropEqualTrue = numeric(),
      stringsAsFactors = FALSE
    )
    for (t_idx in 1:num_T) {
      nf_vec <- results$num_factors[t_idx, ]
      numfactor_table <- rbind(numfactor_table, data.frame(
        T = Tgrid[t_idx],
        MeanNumFactors = mean(nf_vec, na.rm = TRUE),
        SdNumFactors = sd(nf_vec, na.rm = TRUE),
        PropEqualTrue = if (all(is.na(nf_vec))) NA else mean(nf_vec == R_true, na.rm = TRUE)
      ))
    }
    filename_numf <- sprintf("%s_numf_summary_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
    write.csv(numfactor_table, filename_numf, row.names = FALSE)
    cat(sprintf("Estimated-factor summary saved to: %s\n", filename_numf))

    # ============================================================
    # SAVE PER-METHOD ESTIMATES (SEPARATE FILES)
    # ============================================================
    cat("\nSaving per-method estimates...\n")

    build_long_estimates <- function(arr, method_name, se_arr = NULL) {
      df_list <- list()
      idx <- 1
      for (t_idx in 1:num_T) {
        for (j in 1:dimbeta) {
          est_vec <- arr[t_idx, , j]
          se_vec <- NULL
          if (!is.null(se_arr)) {
            se_vec <- se_arr[t_idx, , j]
          }
          df_list[[idx]] <- data.frame(
            T = rep(Tgrid[t_idx], length(est_vec)),
            Simulation = seq_along(est_vec),
            Covariate = rep(j, length(est_vec)),
            Estimate = as.numeric(est_vec),
            SE = if (is.null(se_vec)) rep(NA_real_, length(est_vec)) else as.numeric(se_vec),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
      }
      out <- do.call(rbind, df_list)
      out$Method <- method_name
      out
    }

    # MLE
    se_arr <- sqrt(results$Var_beta)
    mle_df <- build_long_estimates(results$beta_mle, "MLE", se_arr = se_arr)
    filename_mle <- sprintf("%s_estimates_MLE_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
    write.csv(mle_df, filename_mle, row.names = FALSE)
    cat(sprintf("MLE estimates saved to: %s\n", filename_mle))

    # Analytical BC
    ana_df <- build_long_estimates(results$beta_bc_ana, "Analytical", se_arr = se_arr)
    filename_ana <- sprintf("%s_estimates_Analytical_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
    write.csv(ana_df, filename_ana, row.names = FALSE)
    cat(sprintf("Analytical BC estimates saved to: %s\n", filename_ana))

    # Jackknife (Split-panel)
    jack_df <- build_long_estimates(results$beta_corr, "SplitPJ", se_arr = se_arr)
    filename_jack <- sprintf("%s_estimates_SplitPJ_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
    write.csv(jack_df, filename_jack, row.names = FALSE)
    cat(sprintf("Split-panel Jackknife estimates saved to: %s\n", filename_jack))

    # Bootstrap bias-corrected (mean)
    if (num_bootstrap > 0) {
      boot_df <- build_long_estimates(results$beta_bc_boot, "Boot-Mean", se_arr = se_arr)
      filename_boot <- sprintf("%s_estimates_BootMean_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
      write.csv(boot_df, filename_boot, row.names = FALSE)
      cat(sprintf("Bootstrap (mean) estimates saved to: %s\n", filename_boot))

      # Bootstrap bias-corrected (median)
      boot_med_df <- build_long_estimates(results$beta_bc_boot_median, "Boot-Median", se_arr = se_arr)
      filename_boot_med <- sprintf("%s_estimates_BootMedian_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
      write.csv(boot_med_df, filename_boot_med, row.names = FALSE)
      cat(sprintf("Bootstrap (median) estimates saved to: %s\n", filename_boot_med))

      # Bootstrap raw estimates (each row = one simulation, columns = bootstrap reps)
      boot_all_list <- list()
      idx <- 1
      for (t_idx in 1:num_T) {
        for (j in 1:dimbeta) {
          boot_mat <- results$beta_boot_all[t_idx, , j, ]
          if (length(boot_mat) > 0) {
            boot_df <- as.data.frame(boot_mat)
            colnames(boot_df) <- paste0("b", seq_len(ncol(boot_df)))
            boot_df <- cbind(
              T = rep(Tgrid[t_idx], nrow(boot_df)),
              Simulation = seq_len(nrow(boot_df)),
              Covariate = rep(j, nrow(boot_df)),
              boot_df
            )
            boot_all_list[[idx]] <- boot_df
            idx <- idx + 1
          }
        }
      }
      if (length(boot_all_list) > 0) {
        boot_all_df <- do.call(rbind, boot_all_list)
        filename_boot_all <- sprintf("%s_estimates_BootstrapAll_N%d_S%d_scenario%d_delta%.2f.csv", model_type, N, S, dgp_scenario, delta)
        write.csv(boot_all_df, filename_boot_all, row.names = FALSE)
        cat(sprintf("Bootstrap raw estimates saved to: %s\n", filename_boot_all))
      }
    }

    cat("\nDone!\n")
  }
}

