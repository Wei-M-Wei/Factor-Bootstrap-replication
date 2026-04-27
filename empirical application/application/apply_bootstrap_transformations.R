library(bestNormalize)

# Apply CI transformations to a bootstrap draw matrix.
# Each row is one parameter, each column is one bootstrap replication.

compute_ci_for_row <- function(theta_star,
                               theta_hat,
                               se_hat = NA_real_,
                               alpha = 0.05,
                               min_boot = 10) {
  theta_star <- as.numeric(theta_star)
  theta_star <- theta_star[is.finite(theta_star)]

  out <- list(
    normal_lower = NA_real_,
    normal_upper = NA_real_,
    basic_lower = NA_real_,
    basic_upper = NA_real_,
    log_basic_lower = NA_real_,
    log_basic_upper = NA_real_,
    boxcox_basic_lower = NA_real_,
    boxcox_basic_upper = NA_real_,
    yj_basic_lower = NA_real_,
    yj_basic_upper = NA_real_,
    lambda_boxcox = NA_real_,
    lambda_yj = NA_real_,
    n_boot = length(theta_star)
  )

  if (length(theta_star) < min_boot || !is.finite(theta_hat)) {
    return(as.data.frame(out, check.names = FALSE))
  }

  zcrit <- qnorm(1 - alpha / 2)

  if (is.finite(se_hat) && se_hat > 0) {
    out$normal_lower <- theta_hat - zcrit * se_hat
    out$normal_upper <- theta_hat + zcrit * se_hat
  }

  z_star <- theta_star - theta_hat
  q_low <- as.numeric(quantile(z_star, probs = alpha / 2, na.rm = TRUE))
  q_high <- as.numeric(quantile(z_star, probs = 1 - alpha / 2, na.rm = TRUE))

  out$basic_lower <- theta_hat - q_high
  out$basic_upper <- theta_hat - q_low

  if (theta_hat > 0 && all(theta_star > 0)) {
    phi_hat_log <- log(theta_hat)
    phi_star_log <- log(theta_star)
    z_log <- phi_star_log - phi_hat_log

    q_low_log <- as.numeric(quantile(z_log, probs = alpha / 2, na.rm = TRUE))
    q_high_log <- as.numeric(quantile(z_log, probs = 1 - alpha / 2, na.rm = TRUE))

    lower_log_phi <- phi_hat_log - q_high_log
    upper_log_phi <- phi_hat_log - q_low_log

    out$log_basic_lower <- exp(lower_log_phi)
    out$log_basic_upper <- exp(upper_log_phi)

    bc_fit <- boxcox(theta_star, standardize = FALSE)
    out$lambda_boxcox <- bc_fit$lambda

    phi_star_bc <- predict(bc_fit)
    phi_hat_bc <- predict(bc_fit, newdata = theta_hat)
    z_bc <- phi_star_bc - phi_hat_bc

    q_low_bc <- as.numeric(quantile(z_bc, probs = alpha / 2, na.rm = TRUE))
    q_high_bc <- as.numeric(quantile(z_bc, probs = 1 - alpha / 2, na.rm = TRUE))

    lower_bc_phi <- phi_hat_bc - q_high_bc
    upper_bc_phi <- phi_hat_bc - q_low_bc

    out$boxcox_basic_lower <- as.numeric(predict(bc_fit, newdata = lower_bc_phi, inverse = TRUE))
    out$boxcox_basic_upper <- as.numeric(predict(bc_fit, newdata = upper_bc_phi, inverse = TRUE))
  }

  yj_input <- c(theta_hat, theta_star)
  yj_fit <- yeojohnson(yj_input, standardize = FALSE)
  out$lambda_yj <- yj_fit$lambda

  phi_all_yj <- predict(yj_fit)
  phi_hat_yj <- phi_all_yj[1]
  phi_star_yj <- phi_all_yj[-1]

  z_yj <- phi_star_yj - phi_hat_yj
  q_low_yj <- as.numeric(quantile(z_yj, probs = alpha / 2, na.rm = TRUE))
  q_high_yj <- as.numeric(quantile(z_yj, probs = 1 - alpha / 2, na.rm = TRUE))

  lower_yj_phi <- phi_hat_yj - q_high_yj
  upper_yj_phi <- phi_hat_yj - q_low_yj

  out$yj_basic_lower <- as.numeric(predict(yj_fit, newdata = lower_yj_phi, inverse = TRUE))
  out$yj_basic_upper <- as.numeric(predict(yj_fit, newdata = upper_yj_phi, inverse = TRUE))

  as.data.frame(out, check.names = FALSE)
}

apply_transformations <- function(bootstrap_file = "bootstrap_application_estimate.csv",
                                  output_file = "bootstrap_application_transformed_ci.csv",
                                  point_estimates = NULL,
                                  se_estimates = NULL,
                                  alpha_levels = c(0.05, 0.01),
                                  min_boot = 10) {
  if (!file.exists(bootstrap_file)) {
    stop("Bootstrap file not found: ", bootstrap_file)
  }

  boot_df <- read.csv(bootstrap_file, check.names = FALSE)
  boot_mat <- as.matrix(boot_df)
  storage.mode(boot_mat) <- "numeric"

  n_param <- nrow(boot_mat)

  if (is.null(point_estimates)) {
    point_estimates <- rowMeans(boot_mat, na.rm = TRUE)
  }
  if (is.null(se_estimates)) {
    se_estimates <- apply(boot_mat, 1, sd, na.rm = TRUE)
  }

  if (length(point_estimates) != n_param) {
    stop("Length of point_estimates must match number of rows in bootstrap matrix.")
  }
  if (length(se_estimates) != n_param) {
    stop("Length of se_estimates must match number of rows in bootstrap matrix.")
  }

  if (length(alpha_levels) == 0 || any(!is.finite(alpha_levels)) ||
      any(alpha_levels <= 0 | alpha_levels >= 1)) {
    stop("alpha_levels must contain values strictly between 0 and 1.")
  }

  level_suffix <- function(alpha) {
    paste0("_", format(alpha * 100, trim = TRUE, scientific = FALSE), "pct")
  }

  ci_cols <- c(
    "normal_lower", "normal_upper", "normal_length",
    "basic_lower", "basic_upper", "basic_length",
    "log_basic_lower", "log_basic_upper", "log_basic_length",
    "boxcox_basic_lower", "boxcox_basic_upper", "boxcox_basic_length",
    "yj_basic_lower", "yj_basic_upper", "yj_basic_length"
  )

  res_list <- vector("list", n_param)

  for (i in seq_len(n_param)) {
    parameter_name <- if (!is.null(rownames(boot_mat))) rownames(boot_mat)[i] else paste0("param_", i)

    row_out <- data.frame(
      parameter = parameter_name,
      theta_hat = point_estimates[i],
      se_hat = se_estimates[i],
      n_boot = NA_real_,
      stringsAsFactors = FALSE
    )

    lambda_boxcox <- NA_real_
    lambda_yj <- NA_real_

    for (j in seq_along(alpha_levels)) {
      alpha <- alpha_levels[j]
      row_res <- compute_ci_for_row(
        theta_star = boot_mat[i, ],
        theta_hat = point_estimates[i],
        se_hat = se_estimates[i],
        alpha = alpha,
        min_boot = min_boot
      )

      if (j == 1) {
        row_out$n_boot <- row_res$n_boot
        lambda_boxcox <- row_res$lambda_boxcox
        lambda_yj <- row_res$lambda_yj
      }

      row_res$normal_length <- row_res$normal_upper - row_res$normal_lower
      row_res$basic_length <- row_res$basic_upper - row_res$basic_lower
      row_res$log_basic_length <- row_res$log_basic_upper - row_res$log_basic_lower
      row_res$boxcox_basic_length <- row_res$boxcox_basic_upper - row_res$boxcox_basic_lower
      row_res$yj_basic_length <- row_res$yj_basic_upper - row_res$yj_basic_lower

      suff <- level_suffix(alpha)
      ci_block <- row_res[, ci_cols, drop = FALSE]
      names(ci_block) <- paste0(names(ci_block), suff)
      row_out <- cbind(row_out, ci_block)
    }

    row_out$lambda_boxcox <- lambda_boxcox
    row_out$lambda_yj <- lambda_yj

    res_list[[i]] <- row_out
  }

  out <- do.call(rbind, res_list)

  write.csv(out, output_file, row.names = FALSE)

  cat("Saved transformed CI table to:", normalizePath(output_file), "\n")
  return(out)
}

# Run directly for the current application bootstrap file.
res <- apply_transformations(
  bootstrap_file = "bootstrap_application_estimate.csv",
  output_file = "bootstrap_application_transformed_ci_mle.csv",
  point_estimates = c(0.13919130, 0.06472665, 0.26884910, 0.44014596, -0.75907579),
  alpha_levels = c(0.05, 0.01),
  min_boot = 10
)

print(res)
