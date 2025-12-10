#!/usr/bin/env Rscript

#==============================================================================
# EXAMPLE: PRE-SIMULATION PARAMETER VALIDATION
#
# This script demonstrates how to validate parameter combinations BEFORE
# running an expensive Monte Carlo simulation.
#
# Usage: Rscript example_parameter_validation.R
#==============================================================================

library(tidyverse)
library(MASS)
library(lme4)

# Source the functions
source("pm_functions.R")

cat("\n")
cat("╔" %+% strrep("═", 78) %+% "╗\n")
cat("║" %+% strrep(" ", 78) %+% "║\n")
cat("║  PARAMETER VALIDATION EXAMPLE" %+% strrep(" ", 48) %+% "║\n")
cat("║  Pre-simulation check of all parameter combinations" %+% strrep(" ", 25) %+% "║\n")
cat("║" %+% strrep(" ", 78) %+% "║\n")
cat("╚" %+% strrep("═", 78) %+% "╝\n\n")

#==============================================================================
# STEP 1: DEFINE TRIAL DESIGN (Hybrid)
#==============================================================================

cat("STEP 1: Creating trial design (Hybrid design, 8 measurement weeks)...\n\n")

measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)
n_participants <- 70

trial_design <- expand_grid(
  participant_id = 1:n_participants,
  week = measurement_weeks
) %>%
  mutate(
    path = rep(rep(1:4, length.out = n_participants), length(measurement_weeks)),
    treatment = case_when(
      week %in% c(4, 8) ~ 1,
      week %in% c(9, 10) ~ 1,
      week %in% c(11, 12) ~ 0,
      week %in% c(16, 20) ~ 0,
      TRUE ~ 1
    ),
    expectancy = if_else(week <= 16, 1.0, 0.5),
    timepoint_idx = rep(1:length(measurement_weeks), each = n_participants),
    timepoint_name = paste0("t", week)
  )

cat(sprintf("  • Trial design created: %d participants × %d timepoints = %d rows\n",
            n_participants, length(measurement_weeks), nrow(trial_design)))

#==============================================================================
# STEP 2: DEFINE FIXED MODEL PARAMETERS (Hendrickson values)
#==============================================================================

cat("\nSTEP 2: Setting fixed correlation parameters (Hendrickson et al.)...\n\n")

model_params <- list(
  N = n_participants,
  c.tv = 0.8,          # Time-variant autocorrelation (fixed)
  c.pb = 0.8,          # Pharm biomarker autocorrelation (fixed)
  c.br = 0.8,          # Bio response autocorrelation (fixed)
  c.cf1t = 0.2,        # Same-time cross-correlation (fixed)
  c.cfct = 0.1,        # Different-time cross-correlation (fixed)
  c.bm = 0.3,          # Biomarker moderation (will be overridden)
  carryover_t1half = 1.0  # Half-life for carryover
)

cat("  Fixed parameters:\n")
cat(sprintf("    c.tv   = %.2f (within-factor autocorr)\n", model_params$c.tv))
cat(sprintf("    c.pb   = %.2f (within-factor autocorr)\n", model_params$c.pb))
cat(sprintf("    c.br   = %.2f (within-factor autocorr)\n", model_params$c.br))
cat(sprintf("    c.cf1t = %.2f (same-time cross-corr)\n", model_params$c.cf1t))
cat(sprintf("    c.cfct = %.2f (diff-time cross-corr)\n\n", model_params$c.cfct))

#==============================================================================
# STEP 3: DEFINE RESPONSE PARAMETERS
#==============================================================================

cat("STEP 3: Setting response and baseline parameters...\n\n")

resp_param <- tibble(
  cat = c("time_variant", "pharm_biomarker", "bio_response"),
  max = c(1.0, 1.0, 5.0),
  disp = c(2.0, 2.0, 2.0),
  rate = c(0.3, 0.3, 0.3),
  sd = c(2.0, 2.0, 2.0)
)

baseline_param <- tibble(
  cat = c("biomarker", "baseline"),
  m = c(5.0, 10.0),
  sd = c(2.0, 2.0)
)

cat("  Response parameters set (3 factors)\n")
cat("  Baseline parameters set (biomarker + baseline)\n\n")

#==============================================================================
# STEP 4: CREATE PARAMETER GRID FOR VALIDATION
#==============================================================================

cat("STEP 4: Creating parameter grid to validate...\n\n")

# This is where you define all combinations you want to test
param_grid <- expand_grid(
  n_participants = c(70),
  biomarker_correlation = c(0.0, 0.2, 0.3, 0.4, 0.5, 0.6),  # The swept parameter
  carryover_t1half = c(0, 1.0, 2.0)
)

cat(sprintf("  Parameter grid created: %d combinations\n", nrow(param_grid)))
cat("\n")
cat("  Grid contents:\n")
print(param_grid)
cat("\n")

#==============================================================================
# STEP 5: RUN VALIDATION
#==============================================================================

cat("\nSTEP 5: Running pre-simulation validation...\n")
cat("  This will test each parameter combination for positive definiteness\n")
cat("  and compute condition numbers (numerical stability).\n\n")

validation_result <- validate_parameter_grid(
  param_grid = param_grid,
  trial_design = trial_design,
  model_params = model_params,
  resp_param = resp_param,
  baseline_param = baseline_param,
  verbose = TRUE
)

#==============================================================================
# STEP 6: DETAILED REPORT
#==============================================================================

report_parameter_validation(validation_result, param_grid)

#==============================================================================
# STEP 7: FILTER TO VALID COMBINATIONS (if needed)
#==============================================================================

if (validation_result$n_invalid > 0) {
  cat("STEP 6: Filtering to valid combinations only\n")
  cat(sprintf("  Original grid: %d combinations\n", nrow(param_grid)))
  cat(sprintf("  Valid grid:    %d combinations\n", nrow(validation_result$valid_combinations)))
  cat(sprintf("  Excluded:      %d combinations\n\n", validation_result$n_invalid))

  # Update param_grid to only valid combinations
  param_grid_valid <- validation_result$valid_combinations

  cat("Updated parameter grid:\n")
  print(param_grid_valid)
  cat("\n")

  cat("PROCEED WITH SIMULATION using param_grid_valid\n\n")

} else {
  cat("STEP 6: All combinations are valid!\n")
  cat("  You can proceed directly to simulation with the original param_grid\n\n")

  param_grid_valid <- param_grid
}

#==============================================================================
# STEP 8: SUMMARY AND RECOMMENDATIONS
#==============================================================================

cat("="*80, "\n")
cat("VALIDATION SUMMARY & NEXT STEPS\n")
cat("="*80, "\n\n")

cat(sprintf("✓ Valid combinations ready for simulation: %d\n", validation_result$n_valid))

if (validation_result$n_invalid > 0) {
  cat(sprintf("✗ Invalid combinations excluded:           %d\n", validation_result$n_invalid))
  cat("\n  REASON FOR FAILURES:\n")
  cat("  • biomarker_correlation values were too high\n")
  cat("  • Creating non-positive definite covariance matrices\n")
  cat("\n  SOLUTION:\n")
  cat("  • Reduce biomarker_correlation maximum value\n")
  cat("  • Or reduce c.cf1t / c.cfct correlation parameters\n")
}

cat("\n  NUMERICAL CONDITIONING:\n")
if (length(validation_result$condition_numbers) > 0) {
  cond <- validation_result$condition_numbers
  cat(sprintf("  • Condition number range: [%.0f, %.0f]\n", min(cond), max(cond)))
  cat(sprintf("  • All matrices are %s\n",
              if(max(cond) < 100) "well-conditioned ✓" else "some ill-conditioned ⚠"))
}

cat("\n  READY TO PROCEED?\n")
if (validation_result$n_invalid == 0 && max(validation_result$condition_numbers) < 100) {
  cat("  ✓ YES - All combinations pass validation and are numerically stable\n")
  cat("  → Proceed with: simulation_clustered.R or your simulation script\n")
} else if (validation_result$n_invalid == 0) {
  cat("  ⚠ PROCEED WITH CAUTION - Some matrices are ill-conditioned\n")
  cat("  → May experience numerical instability during simulation\n")
} else {
  cat("  ✗ NO - Adjust parameter grid and re-validate\n")
  cat("  → Use filtered param_grid: param_grid_valid\n")
}

cat("\n" %+% strrep("="*80, 1) %+% "\n\n")

cat("Validation complete. Save this output for reference.\n\n")
