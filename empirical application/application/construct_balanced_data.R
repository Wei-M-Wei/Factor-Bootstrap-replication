# Load the dataset
data <- read.csv("rd-data.txt", stringsAsFactors = FALSE)
data$year
# Display basic information
cat("Original dataset:\n")
cat("Number of observations:", nrow(data), "\n")
cat("Number of firms:", length(unique(data$firm)), "\n")
cat("Year range:", min(data$year), "-", max(data$year), "\n")
cat("Number of years:", length(unique(data$year)), "\n")
cat("\n")

# Identify the full time range
min_year <- min(data$year)
max_year <- max(data$year)
all_years <- min_year:max_year
total_years <- length(all_years)

cat("Full year range:", min_year, "-", max_year, "\n")
cat("Total years in range:", total_years, "\n\n")

# Create balanced panel by keeping firms with complete observations across all years
# Count observations per firm
firm_counts <- table(data$firm)
# Find firms with complete data (observations equal to total years)
complete_firms <- names(firm_counts)[firm_counts == total_years]
# Filter to keep only complete firms
balanced_data <- data[data$firm %in% complete_firms, ]
# Sort by firm and year
balanced_data <- balanced_data[order(balanced_data$firm, balanced_data$year), ]

cat("Balanced dataset (firms with complete data across all years):\n")
cat("Number of observations:", nrow(balanced_data), "\n")
cat("Number of firms:", length(unique(balanced_data$firm)), "\n")
cat("Year range:", min(balanced_data$year), "-", max(balanced_data$year), "\n")
cat("\n")

# Check for missing values
cat("Missing data summary:\n")
missing_summary <- sapply(balanced_data, function(x) sum(is.na(x)))
print(missing_summary)
cat("\n")

# Additional statistics
cat("Variable 'lgrd1_dum' (missing indicator) distribution:\n")
table(balanced_data$lgrd1_dum)
cat("\n")

cat("Variable 'pat_any' (patent indicator) distribution:\n")
table(balanced_data$pat_any)
cat("\n")

# Export balanced dataset
write.csv(balanced_data, "rd-data_balanced.csv", row.names = FALSE)
cat("Balanced dataset exported to: rd-data_balanced.csv\n")

# Display first few rows
cat("\nFirst few rows of balanced dataset:\n")
print(head(balanced_data, 15))

# Summary statistics
cat("\n\nSummary statistics of balanced dataset:\n")
print(summary(balanced_data))
