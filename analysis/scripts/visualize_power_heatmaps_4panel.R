# Visualization: 4-Panel Power Heatmaps
# Each panel: carryover (rows) vs effect size (columns)
# 4 panels: 2 designs × 2 model specifications

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# Load results
load("../output/full_pmsim_analysis_hyb_versus_co.RData")

# Use simulation_summary which contains the aggregated results
results <- simulation_summary

cat("Loaded results:\n")
cat(sprintf("  Total rows: %d\n", nrow(results)))
cat(sprintf("  Designs: %s\n", paste(unique(results$design), collapse=", ")))
cat(sprintf("  Effect sizes: %s\n", paste(unique(results$biomarker_correlation), collapse=", ")))
cat(sprintf("  Carryover t1/2: %s\n", paste(unique(results$carryover_t1half), collapse=", ")))

# Create labeled versions for plotting
plot_data <- results %>%
  mutate(
    design_label = case_when(
      design == "hybrid" ~ "Hybrid (4-path)",
      design == "crossover" ~ "Crossover (2-seq)",
      TRUE ~ design
    ),
    model_label = case_when(
      model_carryover ~ "WITH carryover model",
      !model_carryover ~ "WITHOUT carryover model",
      TRUE ~ "Unknown"
    ),
    # Create factor levels for proper ordering
    effect_size_label = factor(
      sprintf("c.bm = %.2f", biomarker_correlation),
      levels = sprintf("c.bm = %.2f", sort(unique(biomarker_correlation)))
    ),
    carryover_label = factor(
      sprintf("t½ = %.1f wk", carryover_t1half),
      levels = sprintf("t½ = %.1f wk", sort(unique(carryover_t1half)))
    )
  )

# Create 4 separate heatmaps
heatmaps <- list()

for (design_type in c("Hybrid (4-path)", "Crossover (2-seq)")) {
  for (model_type in c("WITHOUT carryover model", "WITH carryover model")) {
    
    # Filter data for this panel
    panel_data <- plot_data %>%
      filter(design_label == design_type, model_label == model_type)
    
    # Create heatmap
    p <- ggplot(panel_data, aes(x = effect_size_label, y = carryover_label, fill = power)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), 
                color = "white", size = 4, fontface = "bold") +
      scale_fill_gradient2(
        low = "#2166ac", 
        mid = "#f7f7f7", 
        high = "#b2182b",
        midpoint = 0.5,
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = scales::percent_format(accuracy = 1),
        name = "Power"
      ) +
      labs(
        title = paste(design_type, "|", model_type),
        x = "Effect Size (biomarker correlation)",
        y = "Carryover Half-Life"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size = 9),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    # Store with informative name
    panel_name <- paste(
      gsub("[^A-Za-z0-9]", "_", design_type),
      gsub("[^A-Za-z0-9]", "_", model_type),
      sep = "_"
    )
    heatmaps[[panel_name]] <- p
  }
}

# Combine into 2×2 grid
combined_plot <- (heatmaps[[1]] | heatmaps[[2]]) / 
                 (heatmaps[[3]] | heatmaps[[4]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Power Analysis: N-of-1 Trial Designs with Carryover Effects",
    subtitle = sprintf("N = %d participants | %d iterations per condition | α = 0.05",
                      unique(results$n_participants)[1],
                      unique(results$n_iterations)[1]),
    caption = "Note: Power shown as percentage. Hybrid design uses 4-path randomization; Crossover uses 2-sequence design.\nCarryover modeled when specified; otherwise unmodeled carryover reduces power.",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

# Save plot
output_file_pdf <- "../output/power_heatmaps_4panel.pdf"
output_file_png <- "../output/power_heatmaps_4panel.png"

ggsave(output_file_pdf, combined_plot, width = 12, height = 10, units = "in")
ggsave(output_file_png, combined_plot, width = 12, height = 10, units = "in", dpi = 300)

cat("\n")
cat(strrep("=", 70), "\n")
cat("VISUALIZATION COMPLETE\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("Saved to:\n"))
cat(sprintf("  %s\n", output_file_pdf))
cat(sprintf("  %s\n", output_file_png))
cat("\n")

# Print summary statistics
cat("Power Summary by Panel:\n\n")
summary_stats <- plot_data %>%
  group_by(design_label, model_label) %>%
  summarize(
    min_power = min(power),
    mean_power = mean(power),
    max_power = max(power),
    .groups = 'drop'
  ) %>%
  mutate(across(c(min_power, mean_power, max_power), ~sprintf("%.1f%%", . * 100)))

print(summary_stats, n = 20)

cat("\n\nType I Error Rates (c.bm = 0):\n\n")
type1_error <- plot_data %>%
  filter(biomarker_correlation == 0) %>%
  select(design_label, model_label, carryover_label, power) %>%
  arrange(design_label, model_label, carryover_label) %>%
  mutate(power = sprintf("%.1f%%", power * 100))

print(type1_error, n = 20)

cat("\n")
