# Reproducibility Check Script
# Run this script to verify that the analysis can be fully reproduced

# Load required packages
library(here)
library(sessioninfo)
library(renv)
library(digest)

# 1. Environment Check ====
cat("=== REPRODUCIBILITY CHECK ===\\n")
cat("Date:", as.character(Sys.Date()), "\\n")
cat("Time:", as.character(Sys.time()), "\\n\\n")

# Check R version
cat("R Version:", R.version.string, "\\n")

# Check package versions
cat("\\n=== PACKAGE ENVIRONMENT ===\\n")
if (file.exists("renv.lock")) {
  cat("renv.lock found - checking package versions\\n")
  renv::status()
} else {
  cat("WARNING: No renv.lock found - package versions not locked\\n")
}

# 2. File Integrity Check ====
cat("\\n=== FILE INTEGRITY CHECK ===\\n")

# Check for required files
required_files <- c(
  "DESCRIPTION",
  "analysis/paper/paper.Rmd",
  "R/utils.R",
  "scripts"
)

missing_files <- c()
for (file in required_files) {
  if (file.exists(here(file))) {
    cat("✓", file, "exists\\n")
  } else {
    cat("✗", file, "MISSING\\n")
    missing_files <- c(missing_files, file)
  }
}

if (length(missing_files) > 0) {
  cat("\\nERROR: Missing required files:\\n")
  cat(paste("-", missing_files, collapse = "\\n"), "\\n")
  stop("Cannot proceed with missing files")
}

# 3. Data Integrity Check ====
cat("\\n=== DATA INTEGRITY CHECK ===\\n")

# Check for data files
data_dir <- here("data", "raw_data")
if (dir.exists(data_dir)) {
  data_files <- list.files(data_dir, recursive = TRUE)
  cat("Found", length(data_files), "data files\\n")
  
  # Calculate checksums for data files
  if (length(data_files) > 0) {
    cat("\\nData file checksums:\\n")
    for (file in data_files) {
      if (file.size(file.path(data_dir, file)) > 0) {
        checksum <- digest::digest(file.path(data_dir, file), file = TRUE)
        cat("-", file, ":", checksum, "\\n")
      }
    }
  }
} else {
  cat("No raw data directory found\\n")
}

# 4. Script Execution Check ====
cat("\\n=== SCRIPT EXECUTION CHECK ===\\n")

# Test that key scripts can be sourced without error
test_scripts <- list.files(here("scripts"), pattern = "\\\\.R$", full.names = TRUE)

if (length(test_scripts) > 0) {
  for (script in test_scripts) {
    script_name <- basename(script)
    cat("Testing", script_name, "... ")
    
    tryCatch({
      # Test syntax without executing
      parse(script)
      cat("✓ Syntax OK\\n")
    }, error = function(e) {
      cat("✗ Syntax Error:", e$message, "\\n")
    })
  }
} else {
  cat("No R scripts found in scripts directory\\n")
}

# 5. Session Information ====
cat("\\n=== SESSION INFORMATION ===\\n")
session_info()

cat("\\n=== REPRODUCIBILITY CHECK COMPLETE ===\\n")
