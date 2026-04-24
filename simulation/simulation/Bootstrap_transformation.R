library(bestNormalize)

# =========================
# 0. Global settings
# =========================

beta_true  <- 0.5
cov_target <- 1

model_list <- c("logit", "probit")
scenario_list <- c(0, 2, 4, 6)
T_list <- c(20, 30, 40)

base_path <- "your path"

# =========================
# 1. Helper functions
# =========================

escape_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("#", "\\\\#", x)
  x <- gsub("\\$", "\\\\$", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  return(x)
}

# =========================
# 2. Core function:
#    compute one coverage_table
# =========================

compute_coverage_table <- function(model_type,
                                   scenario_id,
                                   T_target,
                                   beta_true = 0.5,
                                   cov_target = 1,
                                   base_path) {
  model_type = 'logit'
  scenario_id= 0
  T_target = 20
  # -------------------------
  # Read files
  # -------------------------
  mc_file <- file.path(
    base_path,
    paste0(model_type, "_estimates_MLE_N30_S1000_scenario", scenario_id, "_delta1.05.csv")
  )

  boot_file <- file.path(
    base_path,
    paste0(model_type, "_estimates_BootstrapAll_N30_S1000_scenario", scenario_id, "_delta1.05.csv")
  )

  if (!file.exists(mc_file)) {
    stop("MLE file not found: ", mc_file)
  }
  if (!file.exists(boot_file)) {
    stop("Bootstrap file not found: ", boot_file)
  }

  mc <- read.csv(mc_file)
  boot <- read.csv(boot_file)

  # -------------------------
  # Subset MLE table
  # -------------------------
  mc_sub <- mc[
    mc$Method == "MLE" &
      mc$T == T_target &
      mc$Covariate == cov_target,
    c("Simulation", "Estimate", "SE")
  ]

  mc_sub <- mc_sub[order(mc_sub$Simulation), ]

  if (nrow(mc_sub) == 0) {
    warning("No MLE rows found for model=", model_type,
            ", scenario=", scenario_id, ", T=", T_target)
    return(NULL)
  }

  # -------------------------
  # Subset bootstrap table
  # -------------------------
  boot_sub <- boot[
    boot$T == T_target &
      boot$Covariate == cov_target,
  ]

  if (nrow(boot_sub) == 0) {
    warning("No bootstrap rows found for model=", model_type,
            ", scenario=", scenario_id, ", T=", T_target)
    return(NULL)
  }

  boot_cols <- grep("^b[0-9]+$", names(boot_sub), value = TRUE)

  if (length(boot_cols) == 0) {
    warning("No bootstrap columns found for model=", model_type,
            ", scenario=", scenario_id, ", T=", T_target)
    return(NULL)
  }

  boot_sub <- boot_sub[, c("Simulation", boot_cols), drop = FALSE]
  boot_sub <- boot_sub[order(boot_sub$Simulation), , drop = FALSE]

  # -------------------------
  # Match simulations
  # -------------------------
  common_sim <- intersect(mc_sub$Simulation, boot_sub$Simulation)

  mc_match <- mc_sub[mc_sub$Simulation %in% common_sim, ]
  boot_match <- boot_sub[boot_sub$Simulation %in% common_sim, ]

  mc_match <- mc_match[order(mc_match$Simulation), ]
  boot_match <- boot_match[order(boot_match$Simulation), ]

  if (nrow(mc_match) == 0) {
    warning("No matched simulations for model=", model_type,
            ", scenario=", scenario_id, ", T=", T_target)
    return(NULL)
  }

  stopifnot(all(mc_match$Simulation == boot_match$Simulation))

  # -------------------------
  # Coverage rate calculation
  # -------------------------
  n_sim <- nrow(mc_match)

  cover_normal <- numeric(n_sim)
  len_normal <- numeric(n_sim)
  miss_lower_normal <- numeric(n_sim)
  miss_upper_normal <- numeric(n_sim)

  cover_basic <- numeric(n_sim)
  len_basic <- numeric(n_sim)
  miss_lower_basic <- numeric(n_sim)
  miss_upper_basic <- numeric(n_sim)

  cover_log <- numeric(n_sim)
  len_log <- numeric(n_sim)
  miss_lower_log <- numeric(n_sim)
  miss_upper_log <- numeric(n_sim)

  cover_boxcox <- numeric(n_sim)
  len_boxcox <- numeric(n_sim)
  miss_lower_boxcox <- numeric(n_sim)
  miss_upper_boxcox <- numeric(n_sim)

  cover_yj <- numeric(n_sim)
  len_yj <- numeric(n_sim)
  miss_lower_yj <- numeric(n_sim)
  miss_upper_yj <- numeric(n_sim)

  lambda_boxcox_vec <- rep(NA_real_, n_sim)
  lambda_yj_vec <- rep(NA_real_, n_sim)

  for (r in 1:n_sim) {

    theta_hat <- mc_match$Estimate[r]
    se_hat    <- mc_match$SE[r]

    theta_star <- as.numeric(boot_match[r, boot_cols, drop = TRUE])
    theta_star <- theta_star[!is.na(theta_star)]

    if (length(theta_star) < 10 || is.na(theta_hat) || is.na(se_hat)) {
      cover_normal[r] <- NA
      len_normal[r] <- NA
      miss_lower_normal[r] <- NA
      miss_upper_normal[r] <- NA

      cover_basic[r] <- NA
      len_basic[r] <- NA
      miss_lower_basic[r] <- NA
      miss_upper_basic[r] <- NA

      cover_log[r] <- NA
      len_log[r] <- NA
      miss_lower_log[r] <- NA
      miss_upper_log[r] <- NA

      cover_boxcox[r] <- NA
      len_boxcox[r] <- NA
      miss_lower_boxcox[r] <- NA
      miss_upper_boxcox[r] <- NA

      cover_yj[r] <- NA
      len_yj[r] <- NA
      miss_lower_yj[r] <- NA
      miss_upper_yj[r] <- NA

      next
    }

    # 1. Normal CI
    lower_norm <- theta_hat - 1.96 * se_hat
    upper_norm <- theta_hat + 1.96 * se_hat

    cover_normal[r] <- (beta_true >= lower_norm) & (beta_true <= upper_norm)
    len_normal[r]   <- upper_norm - lower_norm
    miss_lower_normal[r] <- (beta_true < lower_norm)
    miss_upper_normal[r] <- (beta_true > upper_norm)

    # 2. Basic bootstrap CI
    z_star <- theta_star - theta_hat

    q_low  <- quantile(z_star, 0.025, na.rm = TRUE)
    q_high <- quantile(z_star, 0.975, na.rm = TRUE)

    lower_basic <- theta_hat - q_high
    upper_basic <- theta_hat - q_low

    cover_basic[r] <- (beta_true >= lower_basic) & (beta_true <= upper_basic)
    len_basic[r]   <- upper_basic - lower_basic
    miss_lower_basic[r] <- (beta_true < lower_basic)
    miss_upper_basic[r] <- (beta_true > upper_basic)

    # 3. Log-transform basic bootstrap CI
    if (is.finite(theta_hat) && theta_hat > 0 && all(theta_star > 0, na.rm = TRUE)) {

      phi_hat_log  <- log(theta_hat)
      phi_star_log <- log(theta_star)

      z_log <- phi_star_log - phi_hat_log

      q_low_log  <- quantile(z_log, 0.025, na.rm = TRUE)
      q_high_log <- quantile(z_log, 0.975, na.rm = TRUE)

      lower_log_phi <- phi_hat_log - q_high_log
      upper_log_phi <- phi_hat_log - q_low_log

      lower_log <- exp(lower_log_phi)
      upper_log <- exp(upper_log_phi)

      cover_log[r] <- (beta_true >= lower_log) & (beta_true <= upper_log)
      len_log[r]   <- upper_log - lower_log
      miss_lower_log[r] <- (beta_true < lower_log)
      miss_upper_log[r] <- (beta_true > upper_log)

    } else {
      cover_log[r] <- NA
      len_log[r] <- NA
      miss_lower_log[r] <- NA
      miss_upper_log[r] <- NA
    }

    # 4. Box-Cox basic bootstrap CI
    if (is.finite(theta_hat) && theta_hat > 0 && all(theta_star > 0, na.rm = TRUE)) {

      bc_fit <- boxcox(theta_star, standardize = FALSE)
      lambda_boxcox_vec[r] <- bc_fit$lambda

      phi_star_bc <- predict(bc_fit)
      phi_hat_bc  <- predict(bc_fit, newdata = theta_hat)

      z_bc_tr <- phi_star_bc - phi_hat_bc

      q_low_bc_tr  <- quantile(z_bc_tr, 0.025, na.rm = TRUE)
      q_high_bc_tr <- quantile(z_bc_tr, 0.975, na.rm = TRUE)

      lower_bc_phi <- phi_hat_bc - q_high_bc_tr
      upper_bc_phi <- phi_hat_bc - q_low_bc_tr

      lower_boxcox <- predict(bc_fit, newdata = lower_bc_phi, inverse = TRUE)
      upper_boxcox <- predict(bc_fit, newdata = upper_bc_phi, inverse = TRUE)

      if (is.finite(lower_boxcox) && is.finite(upper_boxcox)) {
        cover_boxcox[r] <- (beta_true >= lower_boxcox) & (beta_true <= upper_boxcox)
        len_boxcox[r]   <- upper_boxcox - lower_boxcox
        miss_lower_boxcox[r] <- (beta_true < lower_boxcox)
        miss_upper_boxcox[r] <- (beta_true > upper_boxcox)
      } else {
        cover_boxcox[r] <- NA
        len_boxcox[r] <- NA
        miss_lower_boxcox[r] <- NA
        miss_upper_boxcox[r] <- NA
      }

    } else {
      cover_boxcox[r] <- NA
      len_boxcox[r] <- NA
      miss_lower_boxcox[r] <- NA
      miss_upper_boxcox[r] <- NA
    }

    # 5. Yeo-Johnson basic bootstrap CI
    yj_input <- c(theta_hat, theta_star)
    yj_fit <- yeojohnson(yj_input, standardize = FALSE)

    lambda_yj_vec[r] <- yj_fit$lambda

    phi_all_yj <- predict(yj_fit)
    phi_hat_yj <- phi_all_yj[1]
    phi_star_yj <- phi_all_yj[-1]

    z_yj <- phi_star_yj - phi_hat_yj

    q_low_yj  <- quantile(z_yj, 0.025, na.rm = TRUE)
    q_high_yj <- quantile(z_yj, 0.975, na.rm = TRUE)

    lower_yj_phi <- phi_hat_yj - q_high_yj
    upper_yj_phi <- phi_hat_yj - q_low_yj

    lower_yj <- predict(yj_fit, newdata = lower_yj_phi, inverse = TRUE)
    upper_yj <- predict(yj_fit, newdata = upper_yj_phi, inverse = TRUE)

    if (is.finite(lower_yj) && is.finite(upper_yj)) {
      cover_yj[r] <- (beta_true >= lower_yj) & (beta_true <= upper_yj)
      len_yj[r]   <- upper_yj - lower_yj
      miss_lower_yj[r] <- (beta_true < lower_yj)
      miss_upper_yj[r] <- (beta_true > upper_yj)
    } else {
      cover_yj[r] <- NA
      len_yj[r] <- NA
      miss_lower_yj[r] <- NA
      miss_upper_yj[r] <- NA
    }
  }

  coverage_table <- data.frame(
    Method = c("Normal",
               "Basic bootstrap",
               "Log-basic bootstrap",
               "Box-Cox-basic bootstrap",
               "Yeo-Johnson-basic bootstrap"),
    CoverageRate = c(mean(cover_normal, na.rm = TRUE),
                     mean(cover_basic, na.rm = TRUE),
                     mean(cover_log, na.rm = TRUE),
                     mean(cover_boxcox, na.rm = TRUE),
                     mean(cover_yj, na.rm = TRUE)),
    AvgLength = c(mean(len_normal, na.rm = TRUE),
                  mean(len_basic, na.rm = TRUE),
                  mean(len_log, na.rm = TRUE),
                  mean(len_boxcox, na.rm = TRUE),
                  mean(len_yj, na.rm = TRUE)),
    LowerMissRate = c(mean(miss_lower_normal, na.rm = TRUE),
                      mean(miss_lower_basic, na.rm = TRUE),
                      mean(miss_lower_log, na.rm = TRUE),
                      mean(miss_lower_boxcox, na.rm = TRUE),
                      mean(miss_lower_yj, na.rm = TRUE)),
    UpperMissRate = c(mean(miss_upper_normal, na.rm = TRUE),
                      mean(miss_upper_basic, na.rm = TRUE),
                      mean(miss_upper_log, na.rm = TRUE),
                      mean(miss_upper_boxcox, na.rm = TRUE),
                      mean(miss_upper_yj, na.rm = TRUE))
  )

  coverage_table$Model <- model_type
  coverage_table$Scenario <- scenario_id
  coverage_table$T <- T_target
  coverage_table$MeanLambdaBoxCox <- mean(lambda_boxcox_vec, na.rm = TRUE)
  coverage_table$MeanLambdaYJ <- mean(lambda_yj_vec, na.rm = TRUE)

  coverage_table <- coverage_table[, c(
    "Model", "Scenario", "T", "Method", "CoverageRate", "AvgLength",
    "LowerMissRate", "UpperMissRate", "MeanLambdaBoxCox", "MeanLambdaYJ"
  )]

  return(coverage_table)
}

# =========================
# 3. LaTeX table generator
#    one model + one scenario -> one table
# =========================

make_latex_table_for_scenario <- function(df, model_type, scenario_id, digits = 4) {

  subdf <- df[df$Model == model_type & df$Scenario == scenario_id, ]

  if (nrow(subdf) == 0) {
    return(paste0("% No results for ", model_type, ", scenario ", scenario_id))
  }

  subdf <- subdf[, c("T", "Method", "CoverageRate", "AvgLength",
                     "LowerMissRate", "UpperMissRate")]

  subdf$Method <- escape_latex(subdf$Method)
  subdf$CoverageRate  <- sprintf(paste0("%.", digits, "f"), subdf$CoverageRate)
  subdf$AvgLength     <- sprintf(paste0("%.", digits, "f"), subdf$AvgLength)
  subdf$LowerMissRate <- sprintf(paste0("%.", digits, "f"), subdf$LowerMissRate)
  subdf$UpperMissRate <- sprintf(paste0("%.", digits, "f"), subdf$UpperMissRate)

  latex_lines <- c()
  latex_lines <- c(latex_lines, "\\begin{table}[!htbp]")
  latex_lines <- c(latex_lines, "\\centering")
  latex_lines <- c(latex_lines,
                   paste0("\\caption{Coverage results for ", model_type,
                          " model, scenario ", scenario_id, "}"))
  latex_lines <- c(latex_lines,
                   paste0("\\label{tab:", model_type, "_scenario", scenario_id, "}"))
  latex_lines <- c(latex_lines, "\\begin{tabular}{c l c c c c}")
  latex_lines <- c(latex_lines, "\\hline")
  latex_lines <- c(latex_lines,
                   "T & Method & CoverageRate & AvgLength & LowerMissRate & UpperMissRate\\\\")
  latex_lines <- c(latex_lines, "\\hline")

  for (i in 1:nrow(subdf)) {
    row_i <- subdf[i, ]
    line_i <- paste0(
      row_i$T, " & ",
      row_i$Method, " & ",
      row_i$CoverageRate, " & ",
      row_i$AvgLength, " & ",
      row_i$LowerMissRate, " & ",
      row_i$UpperMissRate, "\\\\"
    )
    latex_lines <- c(latex_lines, line_i)
  }

  latex_lines <- c(latex_lines, "\\hline")
  latex_lines <- c(latex_lines, "\\end{tabular}")
  latex_lines <- c(latex_lines, "\\end{table}")

  return(paste(latex_lines, collapse = "\n"))
}

# =========================
# 4. Run all models, scenarios, T
# =========================

all_results <- list()
idx <- 1

for (model_type in model_list) {
  for (sc in scenario_list) {
    for (TT in T_list) {
      cat("Running model =", model_type, ", scenario =", sc, ", T =", TT, "...\n")

      res <- tryCatch(
        compute_coverage_table(
          model_type = model_type,
          scenario_id = sc,
          T_target = TT,
          beta_true = beta_true,
          cov_target = cov_target,
          base_path = base_path
        ),
        error = function(e) {
          cat("Error in model =", model_type,
              ", scenario =", sc, ", T =", TT, ":", e$message, "\n")
          return(NULL)
        }
      )

      if (!is.null(res)) {
        all_results[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
}

results_df <- do.call(rbind, all_results)

print(results_df)

# =========================
# 5. Save CSV by model
# =========================

for (model_type in model_list) {
  subdf <- results_df[results_df$Model == model_type, ]

  write.csv(
    subdf,
    file = file.path(base_path,
                     paste0(model_type, "_coverage_results_all_scenarios.csv")),
    row.names = FALSE
  )

  cat("Saved CSV for", model_type, "to:\n",
      file.path(base_path, paste0(model_type, "_coverage_results_all_scenarios.csv")), "\n")
}

# =========================
# 6. Save LaTeX by model
# =========================

for (model_type in model_list) {

  latex_tables <- lapply(scenario_list, function(sc) {
    make_latex_table_for_scenario(results_df, model_type, sc, digits = 4)
  })

  names(latex_tables) <- paste0("scenario", scenario_list)

  for (nm in names(latex_tables)) {
    cat("\n\n")
    cat("% =========================\n")
    cat("% ", model_type, "_", nm, "\n", sep = "")
    cat("% =========================\n")
    cat(latex_tables[[nm]])
    cat("\n\n")
  }

  tex_file <- file.path(base_path, paste0(model_type, "_coverage_tables_by_scenario.tex"))

  full_tex <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage{lmodern}",
    "\\usepackage{booktabs}",
    "\\begin{document}",
    unlist(latex_tables),
    "\\end{document}"
  )

  cat(full_tex, sep = "\n\n", file = tex_file)

  cat("Saved LaTeX tables for", model_type, "to:\n", tex_file, "\n")
}
