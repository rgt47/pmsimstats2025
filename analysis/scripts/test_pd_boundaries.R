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

# Detailed breakdown by design
cat("Detailed results by design:\n")
cat(strrep("-", 70), "\n\n")

results_summary <- results %>%
  group_by(design, c.bm) %>%
  summarize(
    min_autocorr_tested = min(autocorr),
    max_autocorr_tested = max(autocorr),
    max_valid_autocorr = max(autocorr[is_pd], na.rm = TRUE),
    n_valid = sum(is_pd),
    n_total = n(),
    .groups = "drop"
  ) %>%
  mutate(
    max_valid_autocorr = if_else(is.finite(max_valid_autocorr),
                                  max_valid_autocorr,
                                  NA_real_)
  )

print(results_summary, n = Inf)

# Find optimal combinations for each design
cat("\n", strrep("=", 70), "\n")
cat("RECOMMENDED PARAMETER COMBINATIONS\n")
cat(strrep("=", 70), "\n\n")

for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("%s Design:\n", toupper(design_name)))
  cat(strrep("-", 40), "\n")

  design_results <- results %>% filter(design == design_name)

  # Strategy 1: Keep Hendrickson's autocorr = 0.8, find max c.bm
  hendrickson_autocorr <- design_results %>%
    filter(autocorr == 0.8, is_pd)

  if (nrow(hendrickson_autocorr) > 0) {
    max_cbm_at_08 <- max(hendrickson_autocorr$c.bm)
    cat(sprintf("  Strategy 1: Keep autocorr = 0.8\n"))
    cat(sprintf("    → Maximum c.bm = %.2f\n", max_cbm_at_08))
  } else {
    cat("  Strategy 1: Keep autocorr = 0.8\n")
    cat("    → No valid c.bm (all fail)\n")
  }

  # Strategy 2: Keep c.bm = 0.4, find max autocorr
  cbm_04 <- design_results %>%
    filter(c.bm == 0.4, is_pd)

  if (nrow(cbm_04) > 0) {
    max_autocorr_at_04 <- max(cbm_04$autocorr)
    cat(sprintf("  Strategy 2: Keep c.bm = 0.4\n"))
    cat(sprintf("    → Maximum autocorr = %.2f\n", max_autocorr_at_04))
  } else {
    cat("  Strategy 2: Keep c.bm = 0.4\n")
    cat("    → No valid autocorr (all fail)\n")
  }

  # Strategy 3: Balanced reduction
  balanced <- design_results %>%
    filter(is_pd) %>%
    mutate(total_corr = autocorr + c.bm) %>%
    arrange(desc(total_corr)) %>%
    slice(1)

  if (nrow(balanced) > 0) {
    cat(sprintf("  Strategy 3: Balanced (max sum)\n"))
    cat(sprintf("    → autocorr = %.2f, c.bm = %.2f (sum = %.2f)\n",
                balanced$autocorr, balanced$c.bm, balanced$total_corr))
  }

  cat("\n")
}

# Compare designs
cat(strrep("=", 70), "\n")
cat("DESIGN COMPARISON\n")
cat(strrep("=", 70), "\n\n")

comparison <- results %>%
  group_by(design, c.bm) %>%
  summarize(
    max_valid_autocorr = max(autocorr[is_pd], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = design,
    values_from = max_valid_autocorr,
    names_prefix = "max_autocorr_"
  ) %>%
  mutate(
    across(starts_with("max_autocorr"), ~if_else(is.finite(.), ., NA_real_)),
    difference = max_autocorr_hybrid - max_autocorr_crossover
  )

print(comparison, n = Inf)

cat("\nInterpretation:\n")
cat("  - Positive difference: Hybrid allows higher autocorrelation\n")
cat("  - Negative difference: Crossover allows higher autocorrelation\n")
cat("  - NA: At least one design has no valid autocorrelations\n\n")

# Create visualization
cat(strrep("=", 70), "\n")
cat("Creating visualizations...\n")

# Heatmap for each design
for (design_name in c("hybrid", "crossover")) {
  design_data <- results %>% filter(design == design_name)

  p <- ggplot(design_data, aes(x = c.bm, y = autocorr, fill = is_pd)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("FALSE" = "#d73027", "TRUE" = "#1a9850"),
      labels = c("FALSE" = "Non-PD (Invalid)", "TRUE" = "PD (Valid)"),
      name = "Matrix Status"
    ) +
    scale_x_continuous(breaks = seq(0.2, 0.4, by = 0.05)) +
    scale_y_continuous(breaks = seq(0.5, 0.8, by = 0.05)) +
    labs(
      title = sprintf("Positive-Definite Boundary Map: %s Design", toupper(design_name)),
      subtitle = "Fixed: c.cf1t = 0.2, c.cfct = 0.1 (Hendrickson values)",
      x = "Biomarker Correlation (c.bm)",
      y = "Autocorrelation (c.tv = c.pb = c.br)",
      caption = sprintf("Design: %s (%d timepoints) | Green = Valid, Red = Invalid",
                       design_name,
                       if_else(design_name == "hybrid", 8, 4))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      legend.position = "bottom"
    ) +
    # Add Hendrickson's values as reference point
    annotate("point", x = 0.2, y = 0.8, size = 5, shape = 21,
             fill = "yellow", color = "black", stroke = 2) +
    annotate("text", x = 0.2, y = 0.82,
             label = "Hendrickson\n(c.bm=0.2, autocorr=0.8)",
             size = 3, fontface = "bold")

  ggsave(sprintf("../output/pd_boundary_map_%s.png", design_name),
         p, width = 10, height = 8, dpi = 300)
  ggsave(sprintf("../output/pd_boundary_map_%s.pdf", design_name),
         p, width = 10, height = 8)
}

# Combined comparison plot
p_combined <- ggplot(results, aes(x = c.bm, y = autocorr, fill = is_pd)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ design, labeller = labeller(design = c(
    "hybrid" = "Hybrid (8 timepoints)",
    "crossover" = "Crossover (4 timepoints)"
  ))) +
  scale_fill_manual(
    values = c("FALSE" = "#d73027", "TRUE" = "#1a9850"),
    labels = c("FALSE" = "Non-PD (Invalid)", "TRUE" = "PD (Valid)"),
    name = "Matrix Status"
  ) +
  scale_x_continuous(breaks = seq(0.2, 0.4, by = 0.05)) +
  scale_y_continuous(breaks = seq(0.5, 0.8, by = 0.05)) +
  labs(
    title = "Positive-Definite Boundary Comparison",
    subtitle = "Fixed: c.cf1t = 0.2, c.cfct = 0.1 (Hendrickson values)",
    x = "Biomarker Correlation (c.bm)",
    y = "Autocorrelation (c.tv = c.pb = c.br)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("../output/pd_boundary_map_comparison.png",
       p_combined, width = 14, height = 7, dpi = 300)
ggsave("../output/pd_boundary_map_comparison.pdf",
       p_combined, width = 14, height = 7)

cat("✓ Visualizations saved:\n")
cat("  - ../output/pd_boundary_map_hybrid.png/.pdf\n")
cat("  - ../output/pd_boundary_map_crossover.png/.pdf\n")
cat("  - ../output/pd_boundary_map_comparison.png/.pdf\n\n")

# Save results
save(results, results_summary, comparison, test_grid,
     file = "../output/pd_boundary_test_results.RData")

cat("✓ Results saved to: ../output/pd_boundary_test_results.RData\n")

cat("\n", strrep("=", 70), "\n")
cat("✓ TESTING COMPLETE\n")
cat(strrep("=", 70), "\n")
