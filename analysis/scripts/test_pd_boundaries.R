#!/usr/bin/env Rscript
# test_pd_boundaries.R
#
# Systematically test correlation parameter combinations to find
# the boundary of positive-definiteness for BOTH designs

# Suppress package startup messages
suppressPackageStartupMessages({
  library(tidyverse)
})

# Source the functions
source("pm_functions.R")

cat("Testing Positive-Definiteness Boundaries\n")
cat(strrep("=", 70), "\n\n")

# Response parameters (shared) - must be data frame with cat, max, disp, rate, sd
resp_param <- tibble(
  cat = c("time_variant", "pharm_biomarker", "bio_response"),
  max = c(15.0, 15.0, 15.0),
  disp = c(1.0, 1.0, 1.0),
  rate = c(0.1, 0.1, 0.1),
  sd = c(1.5, 1.5, 1.5)  # Standard deviations
)

# Baseline parameters (shared) - must be data frame with cat, sd, and m (mean)
baseline_param <- tibble(
  cat = c("biomarker", "baseline"),
  sd = c(1.0, 1.0),
  m = c(0.0, 0.0)  # Means (set to 0 for baseline)
)

# Create HYBRID trial design (8 timepoints, 4 periods) - SINGLE PARTICIPANT
hybrid_design <- tibble(
  participant_id = 1,
  period = rep(1:4, each = 2),
  week = 1:8,
  timepoint_name = paste0("W", 1:8),
  treatment = c("A", "B", "A", "B", "A", "B", "A", "B"),
  tod = c(1, 0, 1, 0, 1, 0, 1, 0),
  path = 1,
  e = 1.0,  # Expectancy factor (typically 1.0)
  t_wk = 1,  # Time per week
  tsd = 0    # Time since discontinuation (simplified)
) %>%
  mutate(tpb = cumsum(e))  # Cumulative sum of expectancy

# Create CROSSOVER trial design (4 timepoints, 2 periods) - SINGLE PARTICIPANT
crossover_design <- tibble(
  participant_id = 1,
  period = rep(1:2, each = 2),
  week = 1:4,
  timepoint_name = paste0("W", 1:4),
  treatment = c("A", "B", "A", "B"),
  tod = c(1, 0, 1, 0),
  path = 1,
  e = 1.0,  # Expectancy factor (typically 1.0)
  t_wk = 1,  # Time per week
  tsd = 0    # Time since discontinuation (simplified)
) %>%
  mutate(tpb = cumsum(e))  # Cumulative sum of expectancy

# Test grid
test_grid <- expand_grid(
  # Test both designs
  design = c("hybrid", "crossover"),
  # Test autocorrelations from 0.8 down to 0.5
  autocorr = seq(0.8, 0.5, by = -0.05),
  # Test biomarker correlation from 0.4 down to 0.2
  c.bm = seq(0.4, 0.2, by = -0.05),
  # Fixed carryover
  carryover = 0
)

cat(sprintf("Testing %d parameter combinations (both designs)...\n\n", nrow(test_grid)))

# Test each combination
results <- test_grid %>%
  rowwise() %>%
  mutate(
    is_pd = {
      # Select appropriate design
      trial_design <- if (design == "hybrid") hybrid_design else crossover_design

      # Create model params with current values
      test_params <- list(
        c.tv = autocorr,
        c.pb = autocorr,
        c.br = autocorr,
        c.bm = c.bm,
        c.cf1t = 0.2,  # Keep cross-correlations fixed (Hendrickson)
        c.cfct = 0.1,  # Keep cross-correlations fixed (Hendrickson)
        N = 10,
        carryover_t1half = carryover
      )

      # Try to build sigma matrix
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

# Summary by design and biomarker correlation
cat("Results by Design and Biomarker Correlation:\n")
cat(strrep("-", 70), "\n\n")

for (design_name in c("hybrid", "crossover")) {
  cat(sprintf("%s Design:\n", toupper(design_name)))

  design_results <- results %>% filter(design == design_name)

  for (bm_val in unique(design_results$c.bm)) {
    bm_results <- design_results %>% filter(c.bm == bm_val)

    # Find maximum autocorrelation that works
    valid_autocorr <- bm_results %>% filter(is_pd) %>% pull(autocorr)

    if (length(valid_autocorr) > 0) {
      max_valid_autocorr <- max(valid_autocorr)
      cat(sprintf("  c.bm = %.2f: ✓ Works with autocorr <= %.2f\n",
                  bm_val, max_valid_autocorr))
    } else {
      cat(sprintf("  c.bm = %.2f: ✗ No valid autocorrelations found\n",
                  bm_val))
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
