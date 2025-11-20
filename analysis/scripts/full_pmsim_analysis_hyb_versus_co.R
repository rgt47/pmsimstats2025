# Simplified N-of-1 Trial Simulation - Constant Effect Model
# Maximum simplification for learning

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(MASS)
  library(conflicted)
})

# Resolve conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(lmerTest::lmer)

# =============================================================================
# PARAMETERS
# =============================================================================

n_participants <- 70
n_iterations <- 5
treatment_effect <- 5.0        # Constant effect when on treatment
baseline_mean <- 10.0          # Mean baseline response
between_subject_sd <- 3.0      # SD of participant random effects
within_subject_sd <- 2.8       # SD of measurement noise
biomarker_mean <- 5.0          # Mean biomarker value
biomarker_sd <- 2.0            # SD of biomarker

# Parameter grid - what we're testing
param_grid <- expand_grid(
  biomarker_correlation = c(0.3),
  carryover_t1half = c(0, 0.5)
)

# =============================================================================
# CREATE HYBRID DESIGN
# =============================================================================

# Measurement schedule
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

# Assign participants to 4 paths (balanced)
path_assignment <- sample(rep(1:4, length.out = n_participants))

# Create design matrix
hybrid_design <- expand_grid(
  participant_id = 1:n_participants,
  week = measurement_weeks
) %>%
  mutate(
    path = path_assignment[participant_id],
    # Treatment assignment based on path and week
    treatment = case_when(
      week %in% c(4, 8) ~ 1,              # Open-label weeks: all on treatment
      week == 9 ~ 1,                       # Transition week
      week %in% c(10, 11, 12) & path %in% c(1, 2) ~ 1,  # First crossover period
      week %in% c(10, 11, 12) & path %in% c(3, 4) ~ 0,
      week == 16 & path %in% c(1, 3) ~ 1,  # Second crossover period
      week == 16 & path %in% c(2, 4) ~ 0,
      week == 20 & path %in% c(1, 3) ~ 0,  # Third crossover period
      week == 20 & path %in% c(2, 4) ~ 1,
      TRUE ~ NA_real_
    )
  )

# =============================================================================
# RUN SIMULATION
# =============================================================================

cat("Running simulation with constant effect model...\n")
cat("Treatment effect =", treatment_effect, "\n")
cat("Participants =", n_participants, "\n")
cat("Iterations per condition =", n_iterations, "\n\n")

results <- tibble()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf("Condition %d: biomarker=%.1f, carryover=%.1f\n",
              i, params$biomarker_correlation, params$carryover_t1half))

  # Run iterations
  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    # Generate participant-level random effects
    participant_data <- tibble(
      participant_id = 1:n_participants,
      # Biomarker value for each participant
      biomarker = rnorm(n_participants, biomarker_mean, biomarker_sd),
      # Random intercept (between-subject variability)
      random_intercept = rnorm(n_participants, 0, between_subject_sd)
    )

    # Generate trial data
    trial_data <- hybrid_design %>%
      left_join(participant_data, by = "participant_id") %>%
      group_by(participant_id) %>%
      arrange(week) %>%
      mutate(
        # Calculate time since drug stopped (for carryover)
        drug_stopped = treatment == 0 & lag(treatment, default = 1) == 1,
        time_off = cumsum(drug_stopped),

        # Carryover effect: exponential decay when off treatment
        carryover_effect = if (params$carryover_t1half > 0) {
          ifelse(treatment == 0 & time_off > 0,
                 (0.5)^(time_off / params$carryover_t1half),
                 0)
        } else {
          0
        },

        # Treatment effect (constant when on drug)
        # - Full effect when treatment = 1
        # - Partial effect via carryover when treatment = 0
        drug_effect = treatment * treatment_effect +
                      (1 - treatment) * carryover_effect * treatment_effect,

        # Biomarker moderation: biomarker amplifies treatment effect
        # Higher biomarker = stronger response to treatment
        biomarker_moderation = (biomarker - biomarker_mean) / biomarker_sd *
                               params$biomarker_correlation * drug_effect,

        # Generate response
        # Response = baseline + random effect + drug effect + biomarker moderation + noise
        response = baseline_mean +
                   random_intercept +
                   drug_effect +
                   biomarker_moderation +
                   rnorm(n(), 0, within_subject_sd)
      ) %>%
      ungroup()

    # Fit mixed model
    model_result <- tryCatch({
      # Center biomarker for better interpretation
      trial_data <- trial_data %>%
        mutate(bm_centered = biomarker - mean(biomarker))

      if (params$carryover_t1half > 0) {
        model <- lmer(response ~ treatment * bm_centered + week + carryover_effect +
                      (1 | participant_id),
                      data = trial_data)
      } else {
        model <- lmer(response ~ treatment * bm_centered + week +
                      (1 | participant_id),
                      data = trial_data)
      }

      # Extract treatment × biomarker interaction
      coefs <- summary(model)$coefficients
      idx <- which(rownames(coefs) == "treatment:bm_centered")

      if (length(idx) == 0) {
        stop("Interaction term not found")
      }

      # Calculate p-value
      t_val <- coefs[idx, "t value"]
      df_approx <- nrow(trial_data) - nrow(coefs)
      p_value <- 2 * pt(-abs(t_val), df = df_approx)

      tibble(
        iteration = iter,
        biomarker_correlation = params$biomarker_correlation,
        carryover_t1half = params$carryover_t1half,
        effect_size = coefs[idx, "Estimate"],
        se = coefs[idx, "Std. Error"],
        t_value = t_val,
        p_value = p_value,
        significant = p_value < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(iteration = iter, biomarker_correlation = params$biomarker_correlation,
             carryover_t1half = params$carryover_t1half, effect_size = NA,
             se = NA, t_value = NA, p_value = NA, significant = NA)
    })

    results <- bind_rows(results, model_result)
  }
}

# =============================================================================
# SUMMARIZE RESULTS
# =============================================================================

cat("\n", strrep("=", 50), "\n")
cat("RESULTS\n")
cat(strrep("=", 50), "\n\n")

summary_results <- results %>%
  group_by(biomarker_correlation, carryover_t1half) %>%
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(summary_results)

# =============================================================================
# VISUALIZATION
# =============================================================================

library(viridis)

# Adaptive visualization based on data structure
plot_power_results <- function(data) {
  has_carryover <- n_distinct(data$carryover_t1half) > 1
  has_biomarker <- n_distinct(data$biomarker_correlation) > 1

  # Convert to factors
  data <- data %>%
    mutate(
      carryover_t1half = factor(carryover_t1half),
      biomarker_correlation = factor(biomarker_correlation)
    )

  if (has_carryover && has_biomarker) {
    # 2D heatmap
    p <- ggplot(data, aes(x = carryover_t1half, y = biomarker_correlation, fill = power)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), color = "black", size = 5) +
      labs(x = "Carryover Half-life (weeks)", y = "Biomarker Correlation")

  } else if (has_carryover) {
    # Bar chart by carryover
    p <- ggplot(data, aes(x = carryover_t1half, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Carryover Half-life (weeks)", y = "Power")

  } else if (has_biomarker) {
    # Bar chart by biomarker
    p <- ggplot(data, aes(x = biomarker_correlation, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Biomarker Correlation", y = "Power")

  } else {
    # Single bar
    p <- ggplot(data, aes(x = 1, y = power, fill = power)) +
      geom_col(width = 0.4) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 6) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "", y = "Power") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  # Common styling
  p <- p +
    scale_fill_viridis_c(name = "Power", labels = scales::percent, limits = c(0, 1)) +
    labs(title = "Statistical Power: Constant Effect Model",
         subtitle = sprintf("N=%d, treatment effect=%.1f, %d iterations",
                           n_participants, treatment_effect, n_iterations)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  return(p)
}

# Generate and save plot
p <- plot_power_results(summary_results)
print(p)

# Save outputs
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "power_results.pdf"), p, width = 8, height = 6)
save(results, summary_results, file = file.path(output_dir, "simulation_results.RData"))

cat("\nDone! Results saved to", output_dir, "\n")
cat("- power_results.pdf\n")
cat("- simulation_results.RData\n")
