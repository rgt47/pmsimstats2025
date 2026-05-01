# ============================================================================
# PLOT SIMULATION RESULTS
# ============================================================================
# Loads saved simulation results and generates visualization
# Run this to regenerate plots without re-running simulation
#
# Usage: Rscript plot_results.R
# ============================================================================

library(tidyverse)
library(ggplot2)

# Load saved results
output_dir <- "../output"
load(file.path(output_dir, "simulation_clustered_results.RData"))

cat("Loaded simulation results:\n")
cat("  - n_iterations:", n_iterations, "\n")
cat("  - carryover_mode:", carryover_mode, "\n")
cat("  - Total conditions:", nrow(summary_results), "\n\n")

# ============================================================================
# VISUALIZATION - Two Separate Figures
# ============================================================================
# Figure 1: Main designs (Crossover, Hybrid, OL+BDC) - WITH vs WITHOUT carryover
# Figure 2: Open-Label only (separate plot)

# Prepare plot data with labels
plot_data <- summary_results |>
  filter(!is.na(power)) |>
  mutate(
    bm_mod_label = factor(
      sprintf("%.1f", biomarker_moderation),
      levels = rev(sprintf("%.1f", sort(unique(biomarker_moderation))))
    ),
    carryover_label = factor(
      paste0("t\u00bd=", carryover_halflife, " wk"),
      levels = paste0("t\u00bd=", sort(unique(carryover_halflife)), " wk")
    ),
    n_label = factor(
      paste0("N=", n_participants),
      levels = c("N=35", "N=70")
    ),
    design_label = factor(
      case_when(
        design == "crossover" ~ "Crossover",
        design == "hybrid" ~ "Hybrid/N-of-1",
        design == "ol_bdc" ~ "OL+BDC",
        design == "ol" ~ "Open-Label"
      ),
      levels = c("Crossover", "Hybrid/N-of-1", "OL+BDC", "Open-Label")
    ),
    power_pct = round(power * 100),
    ci_lower_pct = round(ci_lower * 100),
    ci_upper_pct = round(ci_upper * 100)
  )

# ==========================================================================
# FIGURE 1: Main designs (excluding Open-Label)
# ==========================================================================
# For t½=0 (no carryover), duplicate results for both WITH and WITHOUT panels
# since the models are equivalent when there's no carryover in the data
main_data_base <- plot_data |>
  filter(design != "ol")

# Create WITHOUT carryover data (includes t½=0 from model_carryover=TRUE)
without_data <- main_data_base |>
  filter(model_carryover == FALSE | (model_carryover == TRUE & carryover_halflife == 0)) |>
  mutate(
    scenario_label = factor(
      paste0("A. WITHOUT Carryover\n(N=", n_participants, ")"),
      levels = c("A. WITHOUT Carryover\n(N=35)", "A. WITHOUT Carryover\n(N=70)")
    )
  )

# Create WITH carryover data
with_data <- main_data_base |>
  filter(model_carryover == TRUE) |>
  mutate(
    scenario_label = factor(
      paste0("B. WITH Carryover\n(N=", n_participants, ")"),
      levels = c("B. WITH Carryover\n(N=35)", "B. WITH Carryover\n(N=70)")
    )
  )

# Combine
main_data <- bind_rows(without_data, with_data) |>
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = c(
        "A. WITHOUT Carryover\n(N=35)", "A. WITHOUT Carryover\n(N=70)",
        "B. WITH Carryover\n(N=35)", "B. WITH Carryover\n(N=70)"
      )
    )
  )

p_main <- ggplot(main_data, aes(x = carryover_label, y = bm_mod_label, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d", power_pct)),
    color = "black", size = 3.0, fontface = "bold",
    vjust = -0.3
  ) +
  geom_text(
    aes(label = sprintf("(%d, %d)", ci_lower_pct, ci_upper_pct)),
    color = "black", size = 1.8,
    vjust = 1.2
  ) +
  scale_fill_gradient2(
    name = "Power (%)",
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1),
    labels = function(x) sprintf("%.0f", x * 100)
  ) +
  facet_grid(
    scenario_label ~ design_label,
    labeller = label_value
  ) +
  labs(
    title = "Statistical Power for Biomarker-Treatment Interaction",
    subtitle = sprintf(
      "N=35 and N=70 | %d iterations | Clustered designs (8 timepoints)",
      n_iterations
    ),
    x = "Carryover Half-Life",
    y = "Biomarker Moderation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 9),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

# ==========================================================================
# FIGURE 2: Open-Label only (separate figure)
# ==========================================================================
ol_data <- plot_data |>
  filter(design == "ol") |>
  mutate(
    scenario_label = factor(
      paste0("N=", n_participants),
      levels = c("N=35", "N=70")
    )
  )

p_ol <- ggplot(ol_data, aes(x = carryover_label, y = bm_mod_label, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d", power_pct)),
    color = "black", size = 3.5, fontface = "bold",
    vjust = -0.3
  ) +
  geom_text(
    aes(label = sprintf("(%d, %d)", ci_lower_pct, ci_upper_pct)),
    color = "black", size = 2.0,
    vjust = 1.2
  ) +
  scale_fill_gradient2(
    name = "Power (%)",
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1),
    labels = function(x) sprintf("%.0f", x * 100)
  ) +
  facet_grid(
    scenario_label ~ .,
    labeller = label_value
  ) +
  labs(
    title = "Statistical Power: Open-Label Design",
    subtitle = sprintf(
      "%d iterations | 8 timepoints | WITH carryover in analysis model",
      n_iterations
    ),
    x = "Carryover Half-Life",
    y = "Biomarker Moderation"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 10),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display plots
print(p_main)
print(p_ol)

# Save outputs
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Save main figure (3 designs)
ggsave(
  file.path(output_dir, sprintf("power_main_%s.pdf", timestamp)),
  p_main,
  width = 12, height = 10
)
ggsave(
  file.path(output_dir, sprintf("power_main_%s.png", timestamp)),
  p_main,
  width = 12, height = 10, dpi = 300
)

# Save Open-Label figure (separate)
ggsave(
  file.path(output_dir, sprintf("power_openlabel_%s.pdf", timestamp)),
  p_ol,
  width = 6, height = 5
)
ggsave(
  file.path(output_dir, sprintf("power_openlabel_%s.png", timestamp)),
  p_ol,
  width = 6, height = 5, dpi = 300
)

cat("\nPlots saved:\n")
cat(sprintf("  - power_main_%s.pdf/png\n", timestamp))
cat(sprintf("  - power_openlabel_%s.pdf/png\n", timestamp))
