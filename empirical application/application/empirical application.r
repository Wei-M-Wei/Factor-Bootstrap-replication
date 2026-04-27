rm(list = ls())
source("estimate function.R")
library(NNR)

# Load the R&D dataset
df <- read.csv("rd-data_balanced_cleaned.csv", stringsAsFactors = FALSE)
cat("Dataset dimensions:", dim(df), "\n")
cat("Number of firms:", length(unique(df$firm)), "\n")
cat("Number of years:", length(unique(df$year)), "\n")
cat("Year range:", min(df$year), "-", max(df$year), "\n")
cat("Variable names:\n")
print(names(df))
cat("\n")

# Check data structure
cat("Unique years:\n")
print(sort(unique(df$year)))
cat("\nUnique firms (first 10):\n")
print(head(sort(unique(df$firm)), 10))

# No year filtering needed - data already spans 1981-2001
# All firms already have complete data (balanced panel)
df_select <- df

# Verify balance
for (i in unique(df_select$firm)[1:3]) {
  cat("Firm", i, "years:", paste(sort(df_select[which(df_select$firm==i), ]$year), collapse=", "), "\n")
}

# For this dataset, we keep all firms (no exclusions needed)
df_final <- df_select
cat("\nFinal dataset dimensions:", dim(df_final), "\n")
cat("Final number of firms:", length(unique(df_final$firm)), "\n\n")

# Select variables for analysis
# Using: pat_any as outcome, and covariates: lgspilltec1, lgspillsic1, lgmalspilltec1, lgmalspillsic1, lgrd1, lsales1, lgrd1_dum
vars <- c('firm', 'year', 'pat_any', 'lgspilltec1', 'lgspillsic1', 'lgmalspilltec1', 'lgmalspillsic1', 'lgrd1', 'lsales1', 'lgrd1_dum')

df_vars <- df_final[, vars, drop = FALSE]

# Setup for interactive fixed effects
library(dplyr)
names(df_vars)
cat("Unique years:", length(unique(df_vars$year)), "\n")

# Create unit and time IDs
df_vars$unit_id <- as.integer(factor(df_vars$firm, levels = unique(df_vars$firm)))
df_vars$time_id <- as.integer(factor(df_vars$year, levels = sort(unique(df_vars$year))))

# Define covariates (exclude firm, year, unit_id, time_id, and outcome pat_any)
covs <- setdiff(names(df_vars), c("firm", "year", "unit_id", "time_id", 'pat_any', 'lgmalspilltec1', 'lgmalspillsic1'))
cat("Covariates included:", paste(covs, collapse=", "), "\n\n")

# Prepare dataframe for estimation
df_ready <- df_vars[, c("unit_id", "time_id", 'pat_any', covs), drop = FALSE]

# Keep original copy before normalization
df_ready_orig <- df_ready

# Normalize continuous numeric covariates (z-score), but exclude dummy variables
num_covs_all <- setdiff(names(df_ready), c("unit_id", "time_id", "pat_any"))
num_covs_all <- num_covs_all[sapply(df_ready[num_covs_all], is.numeric)]

# Identify dummy variables (only values 0 and 1)
is_dummy <- sapply(num_covs_all, function(var) {
  unique_vals <- unique(na.omit(df_ready[[var]]))
  length(unique_vals) <= 2 && all(unique_vals %in% c(0, 1))
})

# Separate continuous and dummy variables
dummy_vars <- names(is_dummy[is_dummy])
continuous_vars <- setdiff(num_covs_all, dummy_vars)
num_covs <- continuous_vars  # Use continuous vars for estimation

cat("Dummy variables (NOT normalized):", paste(dummy_vars, collapse = ", "), "\n")
cat("Continuous variables (to be normalized):", paste(continuous_vars, collapse = ", "), "\n\n")

# Normalize only continuous variables
if(length(continuous_vars) > 0){
  df_ready[continuous_vars] <- as.data.frame(scale(df_ready[continuous_vars]))
  message("Normalized continuous covariates:", paste(continuous_vars, collapse = ", "))
} else {
  message("No continuous covariates found to normalize.")
}

# Add dummy variables back to num_covs for analysis (they stay unscaled)
num_covs <- c(continuous_vars, dummy_vars)

cat("\nAll covariates for estimation:", paste(num_covs, collapse = ", "), "\n")
cat("Covariate count:", length(num_covs), "\n\n")
df_ready <- df_ready[, c("unit_id", "time_id", "pat_any", num_covs), drop = FALSE]

# =============================================================================
# EXTRACT OUTCOME AND COVARIATES INTO MATRIX FORM (N×T)
# =============================================================================
cat("=== CONVERTING TO MATRIX FORM (N×T) ===\n")

panel_N <- length(unique(df_ready$unit_id))
panel_T <- length(unique(df_ready$time_id))

cat("Panel dimensions: N (firms) =", panel_N, ", T (years) =", panel_T, "\n")
cat("Total observations:", panel_N * panel_T, "\n\n")

# Create outcome matrix Y_mat (N rows = firms, T cols = years)
Y_mat <- matrix(NA_real_, nrow = panel_N, ncol = panel_T)
for(r in seq_len(nrow(df_ready))){
  ui <- df_ready$unit_id[r]
  ti <- df_ready$time_id[r]
  Y_mat[ui, ti] <- df_ready$pat_any[r]
}

cat("✓ Outcome matrix Y_mat created: ", nrow(Y_mat), "×", ncol(Y_mat), "\n")
cat("  - Outcome variable: pat_any (patent indicator)\n")
cat("  - No missing values in outcome matrix:", all(!is.na(Y_mat)), "\n\n")

# Create covariate list: each element is an N×T matrix
X_list <- lapply(num_covs, function(varname){
  M <- matrix(NA_real_, nrow = panel_N, ncol = panel_T)
  for(r in seq_len(nrow(df_ready))){
    ui <- df_ready$unit_id[r]
    ti <- df_ready$time_id[r]
    M[ui, ti] <- df_ready[[varname]][r]
  }
  return(M)
})
names(X_list) <- num_covs

cat("✓ Covariate list X_list created with", length(X_list), "variables:\n")
for(i in seq_along(X_list)){
  var_name <- names(X_list)[i]
  is_complete <- all(!is.na(X_list[[i]]))
  cat("  -", var_name, "(", nrow(X_list[[i]]), "×", ncol(X_list[[i]]), ") - Complete:", is_complete, "\n")
}
cat("\n✓ Matrix conversion complete!\n")
cat("  - Use Y_mat for outcome\n")
cat("  - Use X_list[[var_name]] for each covariate\n")
cat("  - Use X_list for all covariates\n\n")


# =============================================================================
# ESTIMATION
# =============================================================================

est = NNRPanel_estimate(data_frame = as.matrix(df_ready), func = "probit", R_max = 3, delta = 1.05, iter_max = 10000, tol = 1e-8, delta_1 = 0.5)
est_ana = est$beta_corr_data
se_ana = est$std_beta_raw_data
# Calculate z-statistics, two-sided p-values and significance stars for est_ana
z_stats_ana <- est_ana / se_ana
z_stats_ana[se_ana == 0] <- NA
p_values_ana <- 2 * pnorm(-abs(z_stats_ana))
sig_ana <- ifelse(is.na(p_values_ana), NA,
        ifelse(p_values_ana < 0.001, '***',
          ifelse(p_values_ana < 0.01, '**',
            ifelse(p_values_ana < 0.05, '*',
              ifelse(p_values_ana < 0.1, '.', '')))))

est_ana_pvals <- data.frame(
  estimate = as.numeric(est_ana),
  se = as.numeric(se_ana),
  z = as.numeric(z_stats_ana),
  p = as.numeric(p_values_ana),
  sig = sig_ana,
  row.names = num_covs,
  stringsAsFactors = FALSE
)

print(est_ana_pvals)

# Bootstrap and additional analysis use the matrices created above (Y_mat, X_list)

# Call bootstrap using panel-shaped objects and correct N/T
boot_res <- tryCatch({
  bootstrap_bias_correction(
    Y = Y_mat,
    X_list = X_list,
    beta_mle = est$beta_fe[1:length(num_covs)],
    lambda_mle = est$L_fe,
    f_mle = est$R_fe,
    model_type = "probit",
    num_bootstrap = 399,
    N = panel_N,
    T = panel_T,
    R = est$num_factor_est,
    use_parallel = TRUE,
    dimbeta = length(num_covs),
    alpha = 0.01
  )
}, error = function(e) NULL)
est_boot = boot_res$beta_bc_median
se_boot = apply(boot_res$beta_boot, 1, sd)

# To get 5% significant level, please run the following and check its CIs.
boot_res <- tryCatch({
  bootstrap_bias_correction(
    Y = Y_mat,
    X_list = X_list,
    beta_mle = est$beta_fe[1:length(num_covs)],
    lambda_mle = est$L_fe,
    f_mle = est$R_fe,
    model_type = "probit",
    num_bootstrap = 399,
    N = panel_N,
    T = panel_T,
    R = est$num_factor_est,
    use_parallel = TRUE,
    dimbeta = length(num_covs),
    alpha = 0.05
  )
}, error = function(e) NULL)
