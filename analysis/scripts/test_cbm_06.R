#!/usr/bin/env Rscript
# test_cbm_06.R
# Test if c.bm = 0.6 works with autocorr = 0.6

suppressPackageStartupMessages({
  library(tidyverse)
})

source("full_pmsim_analysis_hyb_versus_co.R")

cat("Testing c.bm = 0, 0.3, 0.6 with autocorr = 0.6\n")
cat(strrep("=", 70), "\n\n")

# Test all three c.bm values
for (cbm_val in c(0, 0.3, 0.6)) {
  cat(sprintf("Testing c.bm = %.1f:\n", cbm_val))

  for (design_name in c("hybrid", "crossover")) {
    trial_design <- if (design_name == "hybrid") {
      designs$hybrid$design_paths[[1]]
    } else {
      designs$crossover$design_paths[[1]]
    }

    test_params <- model_params
    test_params$c.tv <- 0.6
    test_params$c.pb <- 0.6
    test_params$c.br <- 0.6
    test_params$c.bm <- cbm_val
    test_params$carryover_t1half <- 0

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
      cat(sprintf("  ✓ %s: VALID\n", design_name))
    } else {
      cat(sprintf("  ✗ %s: FAILED\n", design_name))
    }
  }
  cat("\n")
}

cat("Test complete!\n")
