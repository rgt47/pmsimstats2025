#!/usr/bin/env Rscript
# test_boundary_autocorr_06.R
# Find maximum c.bm with autocorr = 0.6

suppressPackageStartupMessages({
  library(tidyverse)
})

source("full_pmsim_analysis_hyb_versus_co.R")

cat("Finding PD boundary with autocorr = 0.6\n")
cat(strrep("=", 70), "\n\n")

# Test c.bm from 0.2 to 0.7 in 0.01 increments
test_grid <- expand_grid(
  design = c("hybrid", "crossover"),
  autocorr = 0.6,
  c.bm = seq(0.2, 0.7, by = 0.01),
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

# Find boundary
cat("BOUNDARY AT AUTOCORR = 0.6:\n")
cat(strrep("-", 70), "\n\n")

for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("%s Design:\n", toupper(design_name)))

  design_results <- results %>%
    filter(design == design_name, autocorr == 0.6)

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
        cat(sprintf("  carryover t1/2 = %.1f: Maximum c.bm = %.2f (all work)\n",
                    carry, max_cbm))
      }
    } else {
      cat(sprintf("  carryover t1/2 = %.1f: All fail\n", carry))
    }
  }
  cat("\n")
}

cat(strrep("=", 70), "\n")
cat("Test complete!\n")
