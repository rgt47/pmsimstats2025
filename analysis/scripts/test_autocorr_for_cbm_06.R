#!/usr/bin/env Rscript
# test_autocorr_for_cbm_06.R
# Find minimum autocorrelation needed for c.bm = 0.6

suppressPackageStartupMessages({
  library(tidyverse)
})

source("full_pmsim_analysis_hyb_versus_co.R")

cat("Finding minimum autocorr for c.bm = 0.6\n")
cat(strrep("=", 70), "\n\n")

# Test autocorr from 0.3 to 0.6 in 0.05 increments
test_grid <- expand_grid(
  design = c("hybrid", "crossover"),
  autocorr = seq(0.3, 0.6, by = 0.05),
  c.bm = 0.6,
  carryover = c(0, 1.0, 2.0)
)

results <- test_grid %>%
  rowwise() %>%
  mutate(
    is_pd = {
      trial_design <- if (design == "hybrid") {
        designs$hybrid$design_paths[[1]]
      } else {
        designs$crossover$design_paths[[1]]
      }

      test_params <- model_params
      test_params$c.tv <- autocorr
      test_params$c.pb <- autocorr
      test_params$c.br <- autocorr
      test_params$c.bm <- c.bm
      test_params$carryover_t1half <- carryover

      sigma_result <- build_sigma_matrix(
        test_params,
        resp_param,
        baseline_param,
        trial_design,
        factor_types = c("time_variant", "pharm_biomarker", "bio_response"),
        factor_abbreviations = c("tv", "pb", "br"),
        verbose = FALSE
      )

      !is.null(sigma_result)
    }
  ) %>%
  ungroup()

# Find minimum valid autocorr for c.bm = 0.6
cat("MINIMUM AUTOCORR FOR c.bm = 0.6:\n")
cat(strrep("-", 70), "\n\n")

for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("%s Design:\n", toupper(design_name)))

  design_results <- results %>%
    filter(design == design_name, c.bm == 0.6)

  for (carry in unique(design_results$carryover)) {
    carry_results <- design_results %>% filter(carryover == carry)

    valid_autocorr <- carry_results %>% filter(is_pd) %>% pull(autocorr)

    if (length(valid_autocorr) > 0) {
      min_autocorr <- min(valid_autocorr)
      max_fail <- carry_results %>% filter(!is_pd) %>% pull(autocorr)
      max_fail_autocorr <- if (length(max_fail) > 0) max(max_fail) else NA

      if (!is.na(max_fail_autocorr)) {
        cat(sprintf("  carryover t1/2 = %.1f: Minimum autocorr = %.2f (fails at %.2f)\n",
                    carry, min_autocorr, max_fail_autocorr))
      } else {
        cat(sprintf("  carryover t1/2 = %.1f: Minimum autocorr = %.2f (all work)\n",
                    carry, min_autocorr))
      }
    } else {
      cat(sprintf("  carryover t1/2 = %.1f: All fail (need autocorr < %.2f)\n",
                  carry, min(carry_results$autocorr)))
    }
  }
  cat("\n")
}

cat(strrep("=", 70), "\n")
cat("\nRECOMMENDATION:\n")
cat("To use c.bm = 0, 0.3, 0.6, set autocorrelations to:\n")

min_valid <- results %>%
  filter(is_pd) %>%
  pull(autocorr) %>%
  min()

if (length(min_valid) > 0) {
  cat(sprintf("  autocorr >= %.2f\n", min_valid))
} else {
  cat("  No valid autocorr found in tested range (0.3 to 0.6)\n")
  cat("  Need to test lower values\n")
}

cat("\nTest complete!\n")
