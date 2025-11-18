#!/usr/bin/env Rscript
# test_new_params.R
# Quick test: autocorr = 0.6, c.bm = 0.6

suppressPackageStartupMessages({
  library(tidyverse)
})

# Source the simulation script
source("full_pmsim_analysis_hyb_versus_co.R")

cat("Testing new parameters:\n")
cat("  - autocorr = 0.6 (reduced from 0.8)\n")
cat("  - c.bm = 0.6 (increased from 0.34)\n\n")

# Test parameters
test_params <- model_params
test_params$c.tv <- 0.6
test_params$c.pb <- 0.6
test_params$c.br <- 0.6
test_params$c.bm <- 0.6
test_params$carryover_t1half <- 0

# Test both designs
for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("Testing %s design:\n", design_name))

  trial_design <- if (design_name == "hybrid") {
    designs$hybrid$design_paths[[1]]
  } else {
    designs$crossover$design_paths[[1]]
  }

  # Test with different carryover values
  for (carry in c(0, 1.0, 2.0)) {
    test_params$carryover_t1half <- carry

    sigma_result <- build_sigma_matrix(
      test_params,
      resp_param,
      baseline_param,
      trial_design,
      factor_types = c("time_variant", "pharm_biomarker", "bio_response"),
      factor_abbreviations = c("tv", "pb", "br"),
      verbose = TRUE
    )

    if (!is.null(sigma_result)) {
      cat(sprintf("  ✓ carryover t1/2 = %.1f: VALID (PD)\n", carry))
    } else {
      cat(sprintf("  ✗ carryover t1/2 = %.1f: FAILED (non-PD)\n", carry))
    }
  }
  cat("\n")
}

cat("Test complete!\n")
