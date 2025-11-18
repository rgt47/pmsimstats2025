#!/usr/bin/env Rscript
# test_pd_boundaries.R
#
# Systematically test correlation parameter combinations to find
# the boundary of positive-definiteness for BOTH designs
# USES EXACT SAME DESIGN STRUCTURE AS ACTUAL SIMULATION

# Suppress package startup messages
suppressPackageStartupMessages({
  library(tidyverse)
})

# Source the simulation script to get the design creation functions
source("full_pmsim_analysis_hyb_versus_co.R")

cat("Testing Positive-Definiteness Boundaries\n")
cat(strrep("=", 70), "\n\n")

# Use the SAME response and baseline parameters as the actual simulation
# These are already defined in the sourced script:
# - resp_param
# - baseline_param
# - model_params

cat("Using parameters from actual simulation:\n")
cat(sprintf("  - resp_param: %s\n", paste(resp_param$cat, collapse = ", ")))
cat(sprintf("  - baseline_param: %s\n", paste(baseline_param$cat, collapse = ", ")))
cat(sprintf("  - Hendrickson correlations: c.tv=%.2f, c.pb=%.2f, c.br=%.2f\n",
            model_params$c.tv, model_params$c.pb, model_params$c.br))
cat(sprintf("  - Cross-correlations: c.cf1t=%.2f, c.cfct=%.2f\n",
            model_params$c.cf1t, model_params$c.cfct))
cat("\n")

# Use the SAME design templates as the actual simulation
# These were already created by the sourced script:
# - designs$hybrid$design_paths[[1]]
# - designs$crossover$design_paths[[1]]

cat("Using design templates from actual simulation:\n")
cat(sprintf("  - Hybrid design: %d timepoints\n",
            nrow(designs$hybrid$design_paths[[1]])))
cat(sprintf("  - Crossover design: %d timepoints\n",
            nrow(designs$crossover$design_paths[[1]])))
cat("\n")

# Test grid - SAME AS ACTUAL SIMULATION with expanded search
test_grid <- expand_grid(
  # Test both designs
  design = c("hybrid", "crossover"),
  # Focus on autocorr = 0.8 (Hendrickson value) plus nearby values
  autocorr = c(0.75, 0.8),
  # Test biomarker correlation from 0.2 to 0.4 in small increments
  c.bm = seq(0.2, 0.4, by = 0.01),
  # Test with and without carryover
  carryover = c(0, 1.0, 2.0)
)

cat(sprintf("Testing %d parameter combinations (both designs)...\n", nrow(test_grid)))
cat("Grid:\n")
cat(sprintf("  - autocorr: %s\n", paste(unique(test_grid$autocorr), collapse = ", ")))
cat(sprintf("  - c.bm: %.2f to %.2f (step %.2f)\n", min(test_grid$c.bm), max(test_grid$c.bm), 0.01))
cat(sprintf("  - carryover: %s\n", paste(unique(test_grid$carryover), collapse = ", ")))
cat("\n")

# Test each combination
results <- test_grid %>%
  rowwise() %>%
  mutate(
    is_pd = {
      # Select appropriate design template (SAME AS ACTUAL SIMULATION)
      trial_design <- if (design == "hybrid") {
        designs$hybrid$design_paths[[1]]
      } else {
        designs$crossover$design_paths[[1]]
      }

      # Create model params with current values (SAME AS ACTUAL SIMULATION)
      test_params <- model_params  # Start with base params
      test_params$c.tv <- autocorr
      test_params$c.pb <- autocorr
      test_params$c.br <- autocorr
      test_params$c.bm <- c.bm
      test_params$carryover_t1half <- carryover
      # c.cf1t and c.cfct already set to Hendrickson values in model_params
      # N already set in model_params

      # Try to build sigma matrix (SAME CALL AS ACTUAL SIMULATION)
      sigma_result <- build_sigma_matrix(
        test_params,
        resp_param,
        baseline_param,
        trial_design,
        factor_types = c("time_variant", "pharm_biomarker", "bio_response"),
        factor_abbreviations = c("tv", "pb", "br"),
        verbose = FALSE
      )

      # Return TRUE if PD, FALSE if NULL (non-PD)
      !is.null(sigma_result)
    }
  ) %>%
  ungroup()

# Find the exact boundary for autocorr = 0.8
cat(strrep("=", 70), "\n")
cat("EXACT BOUNDARY AT AUTOCORR = 0.8 (HENDRICKSON VALUE)\n")
cat(strrep("=", 70), "\n\n")

for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("%s Design:\n", toupper(design_name)))
  cat(strrep("-", 40), "\n")

  design_results <- results %>%
    filter(design == design_name, autocorr == 0.8)

  # Find boundary by carryover level
  for (carry in unique(design_results$carryover)) {
    carry_results <- design_results %>% filter(carryover == carry)

    valid_cbm <- carry_results %>% filter(is_pd) %>% pull(c.bm)

    if (length(valid_cbm) > 0) {
      max_cbm <- max(valid_cbm)
      min_fail <- carry_results %>% filter(!is_pd) %>% pull(c.bm)
      min_fail_cbm <- if (length(min_fail) > 0) min(min_fail) else NA

      if (!is.na(min_fail_cbm)) {
        cat(sprintf("  carryover t1/2 = %.1f: Maximum c.bm = %.2f (fails at %.2f)\n",
                    carry, max_cbm, min_fail_cbm))
      } else {
        cat(sprintf("  carryover t1/2 = %.1f: Maximum c.bm = %.2f (all tested values work)\n",
                    carry, max_cbm))
      }
    } else {
      cat(sprintf("  carryover t1/2 = %.1f: All fail (even c.bm = %.2f)\n",
                  carry, min(carry_results$c.bm)))
    }
  }
  cat("\n")
}

cat(strrep("-", 70), "\n\n")

# Find optimal combinations for each design
cat(strrep("=", 70), "\n")
cat("RECOMMENDED PARAMETER COMBINATIONS\n")
cat(strrep("=", 70), "\n\n")

for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("%s Design:\n", toupper(design_name)))
  cat(strrep("-", 40), "\n")

  design_results <- results %>% filter(design == design_name)

  # Strategy 1: Keep Hendrickson's autocorr = 0.8, find max c.bm
  hendrickson_autocorr <- design_results %>%
    filter(autocorr == 0.8, carryover == 0, is_pd)

  if (nrow(hendrickson_autocorr) > 0) {
    max_cbm_at_08 <- max(hendrickson_autocorr$c.bm)
    cat(sprintf("  Strategy 1: Keep autocorr = 0.8 (no carryover)\n"))
    cat(sprintf("    → Maximum c.bm = %.2f\n", max_cbm_at_08))
  } else {
    cat("  Strategy 1: Keep autocorr = 0.8 (no carryover)\n")
    cat("    → No valid c.bm (all fail)\n")
  }

  # With carryover
  hendrickson_with_carry <- design_results %>%
    filter(autocorr == 0.8, carryover == 1.0, is_pd)

  if (nrow(hendrickson_with_carry) > 0) {
    max_cbm_at_08_carry <- max(hendrickson_with_carry$c.bm)
    cat(sprintf("  Strategy 2: Keep autocorr = 0.8 (with carryover t1/2 = 1.0)\n"))
    cat(sprintf("    → Maximum c.bm = %.2f\n", max_cbm_at_08_carry))
  }

  cat("\n")
}

# Save results
cat(strrep("=", 70), "\n")
save(results, test_grid, file = "../output/pd_boundary_test_results.RData")

cat("✓ Results saved to: ../output/pd_boundary_test_results.RData\n")
cat(strrep("=", 70), "\n")
cat("✓ TESTING COMPLETE\n")
cat(strrep("=", 70), "\n")
