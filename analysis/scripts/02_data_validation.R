# Data Validation Script
# This script performs comprehensive data quality checks

# Load required packages
library(here)
library(dplyr)
library(visdat)
library(naniar)
library(skimr)
library(janitor)
library(palmerpenguins)

# Source utility functions
source(here("R", "utils.R"))

# Load data - use Palmer penguins dataset as example
# Replace with your actual data loading logic
if (file.exists(here("data", "raw_data", "dataset1.csv"))) {
  raw_data <- readr::read_csv(here("data", "raw_data", "dataset1.csv"))
  cat("Loaded custom dataset from data/raw_data/dataset1.csv\\n")
} else {
  # Use Palmer penguins as example dataset
  data(penguins, package = "palmerpenguins")
  raw_data <- penguins
  cat("Using Palmer penguins dataset for validation example\\n")
}

# 1. BASIC DATA STRUCTURE CHECKS ====
cat("=== BASIC DATA STRUCTURE ===\\n")
cat("Dimensions:", dim(raw_data), "\\n")
cat("Variable names:\\n")
print(names(raw_data))

# 2. MISSING DATA ANALYSIS ====
cat("\\n=== MISSING DATA ANALYSIS ===\\n")
# Overall missingness
miss_var_summary(raw_data)

# Missing data patterns
vis_miss(raw_data)

# Missing data heatmap
gg_miss_upset(raw_data)

# 3. DATA TYPE VALIDATION ====
cat("\\n=== DATA TYPE VALIDATION ===\\n")
# Check data types
glimpse(raw_data)

# 4. OUTLIER DETECTION ====
cat("\\n=== OUTLIER DETECTION ===\\n")
# Statistical summary
skim(raw_data)

# 5. CONSISTENCY CHECKS ====
cat("\\n=== CONSISTENCY CHECKS ===\\n")
# Check for duplicate rows
n_duplicates <- nrow(raw_data) - nrow(distinct(raw_data))
cat("Duplicate rows:", n_duplicates, "\\n")

# 6. RANGE VALIDATION ====
cat("\\n=== RANGE VALIDATION ===\\n")
# Penguin-specific range checks
if ("bill_length_mm" %in% names(raw_data)) {
  bill_length_issues <- raw_data$bill_length_mm < 25 | raw_data$bill_length_mm > 70
  cat("Bill length values outside 25-70mm range:", sum(bill_length_issues, na.rm = TRUE), "\\n")
}

if ("body_mass_g" %in% names(raw_data)) {
  mass_issues <- raw_data$body_mass_g < 2000 | raw_data$body_mass_g > 7000
  cat("Body mass values outside 2000-7000g range:", sum(mass_issues, na.rm = TRUE), "\\n")
}

# Species validation
if ("species" %in% names(raw_data)) {
  expected_species <- c("Adelie", "Chinstrap", "Gentoo")
  unexpected_species <- !raw_data$species %in% expected_species
  cat("Unexpected species values:", sum(unexpected_species, na.rm = TRUE), "\\n")
  if (any(unexpected_species, na.rm = TRUE)) {
    cat("Found species:", unique(raw_data$species[unexpected_species]), "\\n")
  }
}

cat("\\n=== DATA VALIDATION COMPLETE ===\\n")
cat("Review the output above for any data quality issues\\n")
