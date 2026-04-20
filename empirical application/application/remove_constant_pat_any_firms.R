# Load balanced dataset
data <- read.csv("rd-data_balanced.csv", stringsAsFactors = FALSE)

cat("=== STEP 1: Load Balanced Dataset ===\n")
cat("Initial observations:", nrow(data), "\n")
cat("Initial firms:", length(unique(data$firm)), "\n\n")

# STEP 2: Calculate variation for pat_any (outcome) per firm
cat("=== STEP 2: Identify Firms with Constant pat_any ===\n")

firm_variation <- data.frame()

for (firm_id in unique(data$firm)) {
  firm_data <- data[data$firm == firm_id, ]

  # Calculate standard deviation for pat_any
  pat_any_sd <- sd(firm_data$pat_any, na.rm = TRUE)

  # Get unique values for pat_any
  unique_values <- unique(firm_data$pat_any)

  # Flag firm if pat_any has zero variation (all same value)
  has_constant_pat_any <- is.na(pat_any_sd) || pat_any_sd == 0

  firm_variation <- rbind(firm_variation, data.frame(
    firm = firm_id,
    pat_any_sd = pat_any_sd,
    has_constant_pat_any = has_constant_pat_any,
    unique_values = paste(sort(unique_values), collapse = ",")
  ))
}

# Count firms with constant pat_any
firms_with_constant_pat_any <- sum(firm_variation$has_constant_pat_any)

cat("Firms with CONSTANT pat_any (all 0s or all 1s):", firms_with_constant_pat_any, "\n\n")

# Show breakdown
constant_pat_any_firms <- firm_variation[firm_variation$has_constant_pat_any, ]
cat("Breakdown:\n")
all_zeros <- sum(constant_pat_any_firms$unique_values == "0")
all_ones <- sum(constant_pat_any_firms$unique_values == "1")
cat("- Firms with pat_any = 0 (always): ", all_zeros, "\n")
cat("- Firms with pat_any = 1 (always): ", all_ones, "\n\n")

# STEP 3: Filter out firms with constant pat_any
firms_to_keep <- firm_variation$firm[!firm_variation$has_constant_pat_any]
clean_data <- data[data$firm %in% firms_to_keep, ]

cat("=== STEP 3: Remove Firms with Constant pat_any ===\n")
cat("Firms removed:", nrow(data) / 21 - nrow(clean_data) / 21, "\n")
cat("Observations removed:", nrow(data) - nrow(clean_data), "\n\n")

# STEP 4: Summary of cleaned dataset
cat("=== FINAL CLEANED DATASET ===\n")
cat("Final observations:", nrow(clean_data), "\n")
cat("Final firms:", length(unique(clean_data$firm)), "\n")
cat("Year range:", min(clean_data$year), "-", max(clean_data$year), "\n\n")

# Verify it's still balanced
observations_per_firm <- table(clean_data$firm)
if (all(observations_per_firm == 21)) {
  cat("✓ Dataset remains BALANCED (each firm has 21 years)\n\n")
} else {
  cat("WARNING: Dataset is NOT balanced\n\n")
}

# Show patent distribution in cleaned data
cat("Patent indicator distribution (cleaned data):\n")
print(table(clean_data$pat_any))
cat("\nFirms with pat_any variation (between 0 and 1):", length(unique(clean_data$firm)), "\n\n")

# Export cleaned dataset
write.csv(clean_data, "rd-data_balanced_cleaned.csv", row.names = FALSE)
cat("Cleaned dataset exported to: rd-data_balanced_cleaned.csv\n\n")

# Show summary statistics of cleaned data
cat("Summary statistics of cleaned dataset:\n")
print(summary(clean_data))

# Show which firms were removed
cat("\n=== FIRMS REMOVED (showing first 15) ===\n")
removed_firms <- firm_variation[firm_variation$has_constant_pat_any, ]
if (nrow(removed_firms) > 0) {
  print_removed <- head(removed_firms, 15)
  for (i in 1:nrow(print_removed)) {
    pat_value <- print_removed$unique_values[i]
    cat("Firm", print_removed$firm[i], "- pat_any constant at:", pat_value, "\n")
  }
  if (nrow(removed_firms) > 15) {
    cat("... and", nrow(removed_firms) - 15, "more firms\n")
  }
}

