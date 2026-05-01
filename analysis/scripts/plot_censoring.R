# ============================================================================
# CENSORING SENSITIVITY ANALYSIS - VISUALIZATION ONLY
# ============================================================================
# Loads saved simulation results and regenerates plots
# Usage: Rscript plot_censoring.R
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Load saved results
output_dir <- "../output"
load(file.path(output_dir, "simulation_censoring_results.RData"))

cat("Loaded simulation results:\n")
cat("  Iterations:", n_iterations, "\n")
cat("  Conditions:", nrow(param_grid), "\n")
cat("  Results rows:", nrow(results), "\n\n")

# ============================================================================
# PREPARE PLOT DATA
# ============================================================================

plot_data <- summary_results |>
  filter(!is.na(power)) |>
  mutate(
    censoring_label = factor(
      censoring_type,
      levels = c("high_dropout", "more_biased", "more_flat", "balanced", "none"),
      labels = c("High dropout", "More biased", "More flat", "Balanced", "No censoring")
    ),
    bm_label = factor(
      biomarker_moderation,
      levels = c(0, 0.3, 0.6)
    ),
    carryover_label = factor(
      paste0("t\u00bd=", carryover_halflife),
      levels = paste0("t\u00bd=", c(0, 0.1, 0.2))
    ),
    design_label = factor(
      case_when(
        design == "crossover" ~ "Crossover",
        design == "hybrid" ~ "N-of-1",
        design == "ol_bdc" ~ "OL+BDC",
        design == "ol" ~ "Open-Label"
      ),
      levels = c("Open-Label", "OL+BDC", "Crossover", "N-of-1")
    ),
    n_label = factor(paste0("N=", n_participants), levels = c("N=35", "N=70")),
    power_pct = round(power * 100)
  )

# ============================================================================
# HENDRICKSON FIGURE 4A STYLE HEATMAP
# ============================================================================
# Layout: Carryover on x-axis, censoring on y-axis
# Faceted by design (columns) and N + biomarker (rows)

p_main <- ggplot(
  plot_data,
  aes(x = carryover_label, y = censoring_label, fill = power)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", power)),
    color = "black", size = 2.2
  ) +
  scale_fill_gradient2(
    name = "Power",
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1)
  ) +
  facet_grid(
    n_label + bm_label ~ design_label,
    labeller = labeller(
      n_label = label_value,
      bm_label = function(x) paste0("BM=", x),
      design_label = label_value
    )
  ) +
  labs(
    title = "Effect of Trial Design, Censoring, and Carryover on Power",
    subtitle = sprintf(
      "%d iterations | Hendrickson-style censoring sensitivity analysis",
      n_iterations
    ),
    x = "Carryover Half-Life",
    y = "Censoring Parameters"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 8),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  )

# ============================================================================
# ALTERNATIVE: SEPARATE PLOTS PER SAMPLE SIZE
# ============================================================================

make_n_plot <- function(data, n_val) {
  plot_subset <- data |>
    filter(n_participants == n_val) |>
    mutate(censoring_label = droplevels(censoring_label))

  ggplot(
    plot_subset,
    aes(x = carryover_label, y = censoring_label, fill = power)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = power_pct),
      color = "black", size = 2.5
    ) +
    scale_fill_gradient2(
      name = "Power (%)",
      low = "#d73027", mid = "#fee08b", high = "#1a9850",
      midpoint = 0.5, limits = c(0, 1),
      labels = scales::percent
    ) +
    scale_y_discrete(drop = FALSE) +
    facet_grid(
      bm_label ~ design_label,
      labeller = labeller(
        bm_label = function(x) paste0("BM=", x),
        design_label = label_value
      )
    ) +
    labs(
      title = sprintf("N = %d", n_val),
      x = "Carryover Half-Life",
      y = "Censoring"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      strip.text = element_text(face = "bold", size = 9),
      strip.text.y = element_text(angle = 0),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      panel.spacing = unit(0.5, "lines")
    )
}

p_n35 <- make_n_plot(plot_data, 35)
p_n70 <- make_n_plot(plot_data, 70)

p_combined <- p_n35 / p_n70 +
  plot_annotation(
    title = "Censoring Sensitivity Analysis: Power by Trial Design",
    subtitle = sprintf(
      "%d iterations | Carryover: 0, 0.1, 0.2 weeks | 5 censoring levels",
      n_iterations
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  )

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

ggsave(
  file.path(output_dir, sprintf("power_censoring_%s.pdf", timestamp)),
  p_main,
  width = 14, height = 16
)
ggsave(
  file.path(output_dir, sprintf("power_censoring_%s.png", timestamp)),
  p_main,
  width = 14, height = 16, dpi = 300
)

ggsave(
  file.path(output_dir, sprintf("power_censoring_combined_%s.pdf", timestamp)),
  p_combined,
  width = 14, height = 14
)
ggsave(
  file.path(output_dir, sprintf("power_censoring_combined_%s.png", timestamp)),
  p_combined,
  width = 14, height = 14, dpi = 300
)

cat("Plots saved to:", output_dir, "\n")
cat(sprintf("  - power_censoring_%s.pdf/png\n", timestamp))
cat(sprintf("  - power_censoring_combined_%s.pdf/png\n", timestamp))
