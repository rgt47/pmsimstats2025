#!/usr/bin/env Rscript
# visualize_hendrickson_style.R
#
# Creates Figure 4-equivalent matching Hendrickson's layout but showing
# BOTH analysis scenarios (WITH vs WITHOUT carryover modeling)
#
# Panel A: No carryover (t1/2 = 0), varying biomarker_correlation
# Panel B: With carryover (varying t1/2), fixed biomarker_correlation = 0.4
#
# Each panel shows:
#   - Hybrid vs Crossover design comparison
#   - WITH carryover model vs WITHOUT carryover model

library(tidyverse)
library(patchwork)

# Load simulation results
load("analysis/output/full_pmsim_analysis_hyb_versus_co.RData")

# ============================================================================
# Panel A: No Carryover, Varying Biomarker Correlation
# ============================================================================

panel_a_data <- simulation_summary %>%
  filter(carryover_t1half == 0) %>%
  mutate(
    design_label = str_to_title(design),
    approach_label = if_else(
      model_carryover,
      "WITH Carryover Model",
      "WITHOUT Carryover Model"
    )
  )

panel_a <- ggplot(
  panel_a_data,
  aes(x = factor(biomarker_correlation), y = design_label, fill = power)
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
  facet_wrap(~ approach_label, ncol = 2) +
  labs(
    title = "A. No Carryover (Half-life = 0)",
    subtitle = "Power by Biomarker Correlation",
    x = "Biomarker Correlation",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = "gray50")
  )

# ============================================================================
# Panel B: With Carryover, Fixed Biomarker Correlation = 0.4
# ============================================================================

panel_b_data <- simulation_summary %>%
  filter(
    carryover_t1half > 0,
    biomarker_correlation == 0.4
  ) %>%
  mutate(
    design_label = str_to_title(design),
    approach_label = if_else(
      model_carryover,
      "WITH Carryover Model",
      "WITHOUT Carryover Model"
    ),
    carryover_label = sprintf("t1/2 = %.1f", carryover_t1half)
  )

panel_b <- ggplot(
  panel_b_data,
  aes(x = carryover_label, y = design_label, fill = power)
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
  facet_wrap(~ approach_label, ncol = 2) +
  labs(
    title = "B. With Carryover (Biomarker Correlation = 0.4)",
    subtitle = "Power by Carryover Half-life",
    x = "Carryover Half-life (weeks)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = "gray50")
  )

# ============================================================================
# Combine Panels
# ============================================================================

combined_figure <- panel_a / panel_b +
  plot_annotation(
    title = "Figure 4-Equivalent: Impact of Carryover Modeling on Statistical Power",
    subtitle = sprintf(
      "Monte Carlo simulation (N=%d participants, %d iterations per combination)",
      unique(simulation_summary$n_participants)[1],
      max(simulation_results$iteration)
    ),
    caption = "Comparison of analysis approaches:\n  • WITH Carryover Model: Controls for carryover as covariate (our enhancement)\n  • WITHOUT Carryover Model: Ignores carryover (Hendrickson approach)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  )

# ============================================================================
# Save Figure
# ============================================================================

ggsave(
  filename = "analysis/output/figure4_equivalent_hendrickson_style.png",
  plot = combined_figure,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = "analysis/output/figure4_equivalent_hendrickson_style.pdf",
  plot = combined_figure,
  width = 14,
  height = 10
)

cat("\n✓ Figure 4-equivalent saved to:\n")
cat("  - analysis/output/figure4_equivalent_hendrickson_style.png\n")
cat("  - analysis/output/figure4_equivalent_hendrickson_style.pdf\n")

# ============================================================================
# Print Summary Statistics
# ============================================================================

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 80) %+% "\n\n")

cat("Panel A (No Carryover, carryover_t1half = 0):\n")
cat(strrep("-", 80) %+% "\n")
print(panel_a_data %>%
  select(design_label, biomarker_correlation, approach_label, power) %>%
  arrange(approach_label, biomarker_correlation, design_label))

cat("\n\nPanel B (With Carryover, biomarker_correlation = 0.4):\n")
cat(strrep("-", 80) %+% "\n")
print(panel_b_data %>%
  select(design_label, carryover_t1half, approach_label, power) %>%
  arrange(approach_label, carryover_t1half, design_label))

# ============================================================================
# Key Findings
# ============================================================================

cat("\n\n" %+% strrep("=", 80) %+% "\n")
cat("KEY FINDINGS\n")
cat(strrep("=", 80) %+% "\n\n")

# Power loss from not modeling carryover
power_loss <- simulation_summary %>%
  filter(carryover_t1half > 0) %>%
  select(design, carryover_t1half, biomarker_correlation, model_carryover, power) %>%
  pivot_wider(
    names_from = model_carryover,
    values_from = power,
    names_prefix = "model_"
  ) %>%
  mutate(power_loss = model_TRUE - model_FALSE) %>%
  arrange(design, biomarker_correlation, carryover_t1half)

cat("Power Lost by NOT Modeling Carryover:\n")
cat(strrep("-", 80) %+% "\n")
print(power_loss %>%
  select(design, biomarker_correlation, carryover_t1half,
         `WITH Model` = model_TRUE, `WITHOUT Model` = model_FALSE, power_loss))

cat("\n\n✓ Visualization complete!\n")
