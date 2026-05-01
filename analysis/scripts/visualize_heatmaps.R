#!/usr/bin/env Rscript
# visualize_heatmaps.R
#
# Creates 4 heatmaps showing power across all parameter combinations:
#   - X-axis: Carryover half-life (0, 1.0, 2.0 weeks)
#   - Y-axis: Biomarker correlation (0.2, 0.34)
#   - Facets: 4 combinations of design × analysis approach
#     1. Hybrid design, WITHOUT carryover model
#     2. Hybrid design, WITH carryover model
#     3. Crossover design, WITHOUT carryover model
#     4. Crossover design, WITH carryover model

# Suppress package startup messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(conflicted)
})

# Resolve conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")

# Load simulation results
load("../output/full_pmsim_analysis_hyb_versus_co.RData")

# Prepare data for heatmaps
heatmap_data <- simulation_summary %>%
  mutate(
    design_label = str_to_title(design),
    approach_label = if_else(
      model_carryover,
      "WITH Carryover Model",
      "WITHOUT Carryover Model"
    ),
    carryover_label = sprintf("%.1f", carryover_t1half),
    cbm_label = sprintf("%.2f", biomarker_correlation)
  )

# Create individual heatmaps for each design × approach combination
create_heatmap <- function(data, design_val, approach_val) {
  plot_data <- data %>%
    filter(design_label == design_val, approach_label == approach_val)

  p <- ggplot(
    plot_data,
    aes(x = carryover_label, y = cbm_label, fill = power)
  ) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", power)), size = 5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee090",
      high = "#1a9850",
      midpoint = 0.80,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      name = "Power"
    ) +
    labs(
      title = sprintf("%s Design\n%s", design_val, approach_val),
      x = "Carryover Half-life (weeks)",
      y = "Biomarker Correlation (c.bm)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right"
    )

  return(p)
}

# Create 4 heatmaps
heatmap1 <- create_heatmap(heatmap_data, "Hybrid", "WITHOUT Carryover Model")
heatmap2 <- create_heatmap(heatmap_data, "Hybrid", "WITH Carryover Model")
heatmap3 <- create_heatmap(heatmap_data, "Crossover", "WITHOUT Carryover Model")
heatmap4 <- create_heatmap(heatmap_data, "Crossover", "WITH Carryover Model")

# Combine into 2x2 grid
combined_heatmaps <- (heatmap1 | heatmap2) / (heatmap3 | heatmap4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Statistical Power Across All Parameter Combinations",
    subtitle = sprintf(
      "Monte Carlo simulation (N=%d participants, %d iterations per combination)",
      unique(simulation_summary$n_participants)[1],
      max(simulation_results$iteration)
    ),
    caption = "Each heatmap shows power to detect biomarker×treatment interaction\nBiomarker correlation levels: 0 (no interaction), 0.3 (moderate), 0.48 (strong)\nFixed parameters: autocorr = 0.6 (reduced from Hendrickson 0.8), treatment effect = 5.0",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  )

# Save figure
ggsave(
  filename = "../output/power_heatmaps_2x2.png",
  plot = combined_heatmaps,
  width = 16,
  height = 12,
  dpi = 300
)

ggsave(
  filename = "../output/power_heatmaps_2x2.pdf",
  plot = combined_heatmaps,
  width = 16,
  height = 12
)

cat("\n✓ Heatmap figures saved to:\n")
cat("  - ../output/power_heatmaps_2x2.png\n")
cat("  - ../output/power_heatmaps_2x2.pdf\n")

# ============================================================================
# Print Summary Statistics
# ============================================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 80), "\n\n", sep = "")

for (design_val in c("Hybrid", "Crossover")) {
  for (approach_val in c("WITHOUT Carryover Model", "WITH Carryover Model")) {
    cat(sprintf("%s Design - %s:\n", design_val, approach_val))
    cat(strrep("-", 80), "\n", sep = "")

    subset_data <- heatmap_data %>%
      filter(design_label == design_val, approach_label == approach_val) %>%
      select(biomarker_correlation, carryover_t1half, power) %>%
      arrange(biomarker_correlation, carryover_t1half)

    print(subset_data)
    cat("\n")
  }
}

# ============================================================================
# Key Findings
# ============================================================================

cat(strrep("=", 80), "\n", sep = "")
cat("KEY FINDINGS\n")
cat(strrep("=", 80), "\n\n", sep = "")

# Compare designs
cat("Design Comparison (averaging across all conditions):\n")
cat(strrep("-", 80), "\n", sep = "")
design_comparison <- simulation_summary %>%
  group_by(design_label = str_to_title(design)) %>%
  summarize(
    mean_power = mean(power),
    min_power = min(power),
    max_power = max(power),
    .groups = "drop"
  )
print(design_comparison)

cat("\n\nAnalysis Approach Comparison (averaging across all conditions):\n")
cat(strrep("-", 80), "\n", sep = "")
approach_comparison <- simulation_summary %>%
  mutate(approach_label = if_else(model_carryover,
                                   "WITH Carryover Model",
                                   "WITHOUT Carryover Model")) %>%
  group_by(approach_label) %>%
  summarize(
    mean_power = mean(power),
    min_power = min(power),
    max_power = max(power),
    .groups = "drop"
  )
print(approach_comparison)

# Power loss from not modeling carryover
cat("\n\nPower Difference (WITH - WITHOUT carryover model):\n")
cat(strrep("-", 80), "\n", sep = "")
power_loss <- simulation_summary %>%
  select(design, carryover_t1half, biomarker_correlation, model_carryover, power) %>%
  pivot_wider(
    names_from = model_carryover,
    values_from = power,
    names_prefix = "model_"
  ) %>%
  mutate(power_diff = model_TRUE - model_FALSE) %>%
  arrange(design, biomarker_correlation, carryover_t1half)

print(power_loss %>%
  select(design, biomarker_correlation, carryover_t1half,
         `WITH Model` = model_TRUE, `WITHOUT Model` = model_FALSE, power_diff))

cat("\n\n✓ Visualization complete!\n")
