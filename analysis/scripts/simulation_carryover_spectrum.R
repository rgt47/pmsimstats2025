#!/usr/bin/env Rscript
# simulation_carryover_spectrum.R
#
# Detailed simulation across fine-grained carryover spectrum
# to quantify the non-monotonic power phenomenon
#
# This script generates empirical evidence for the white paper on
# carryover bias-variance tradeoff

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
})

# Source shared functions
source("pm_functions.R")

set.seed(42)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CARRYOVER SPECTRUM SIMULATION\n")
cat("Quantifying Non-Monotonic Power in OL+BDC Without Carryover Model\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ============================================================================
# Parameters
# ============================================================================

n_iterations <- 100
n_participants <- 70
BR_rate <- 0.5
biomarker_mean <- 5.0
biomarker_sd <- 2.0
true_beta_bm <- 0.3  # True biomarker moderation effect
noise_sd <- 2.0

# Fine-grained carryover half-lives (weeks)
carryover_values <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0)

# OL+BDC measurement schedule
measurement_weeks <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 20)

cat(sprintf("Iterations: %d\n", n_iterations))
cat(sprintf("Participants: %d\n", n_participants))
cat(sprintf("Carryover values: %s\n", paste(carryover_values, collapse = ", ")))
cat("\n")

# ============================================================================
# Data Generation Function
# ============================================================================

generate_olbdc_data <- function(n_participants, weeks, carryover_t1half,
                                 BR_rate, true_beta_bm, noise_sd) {
  # Create participant info
  participant_info <- tibble(
    participant_id = 1:n_participants,
    biomarker = rnorm(n_participants, biomarker_mean, biomarker_sd),
    bm_centered = (biomarker - biomarker_mean) / biomarker_sd
  )

  # Create design matrix
  design <- expand_grid(
    participant_id = 1:n_participants,
    week = weeks
  ) |>
    mutate(
      treatment = if_else(week <= 8, 1, 0),
      expectancy = if_else(week <= 8, 1.0, 0.5)
    ) |>
    left_join(participant_info, by = "participant_id") |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      # Track discontinuation
      just_discontinued = treatment == 0 & lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(just_discontinued, week, NA_real_),
      discontinuation_week = zoo::na.locf(discontinuation_week, na.rm = FALSE),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        week - discontinuation_week,
        0
      ),
      weeks_at_discont = if_else(
        just_discontinued,
        lag(weeks_on_drug, default = 0),
        NA_real_
      ),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE)
    ) |>
    ungroup()

  # Generate true biological response
  design <- design |>
    mutate(
      # Carryover retention factor
      carryover_retention = case_when(
        treatment == 1 ~ 1,
        carryover_t1half == 0 ~ 0,
        TRUE ~ (1/2)^(tsd / carryover_t1half)
      ),
      # True BR with biomarker moderation
      true_BR = case_when(
        treatment == 1 ~ weeks_on_drug * BR_rate * (1 + true_beta_bm * bm_centered),
        carryover_t1half > 0 ~ weeks_at_discont * BR_rate *
          (1 + true_beta_bm * bm_centered) * carryover_retention,
        TRUE ~ 0
      ),
      # Add noise
      noise = rnorm(n(), 0, noise_sd),
      # Observed outcome
      outcome = 50 - true_BR + noise  # Higher baseline, lower = better
    )

  design
}

# ============================================================================
# Analysis Function
# ============================================================================

analyze_data <- function(data, model_carryover, carryover_halflife = 0.2) {
  if (model_carryover) {
    # Analysis WITH carryover modeling
    data <- data |>
      mutate(
        effective_weeks = case_when(
          treatment == 1 ~ weeks_on_drug,
          carryover_halflife > 0 ~ weeks_at_discont * (1/2)^(tsd / carryover_halflife),
          TRUE ~ 0
        )
      )
    predictor <- "effective_weeks"
  } else {
    # Analysis WITHOUT carryover modeling (simple predictor)
    predictor <- "weeks_on_drug"
  }

  # Build formula
  formula_str <- sprintf(
    "outcome ~ bm_centered + %s + bm_centered:%s + (1|participant_id)",
    predictor, predictor
  )

  # Fit model
  model <- tryCatch({
    lmer(as.formula(formula_str), data = data, REML = FALSE)
  }, error = function(e) NULL)

  if (is.null(model)) {
    return(tibble(beta = NA, se = NA, p_value = NA))
  }

  # Extract interaction coefficient
  coef_summary <- summary(model)$coefficients
  interaction_name <- paste0("bm_centered:", predictor)

  if (interaction_name %in% rownames(coef_summary)) {
    beta <- coef_summary[interaction_name, 1]
    se <- coef_summary[interaction_name, 2]
    # Approximate p-value (Satterthwaite would be more accurate)
    t_val <- beta / se
    p_value <- 2 * pt(-abs(t_val), df = n_participants - 4)
  } else {
    beta <- NA
    se <- NA
    p_value <- NA
  }

  tibble(beta = beta, se = se, p_value = p_value)
}

# ============================================================================
# Run Simulation
# ============================================================================

results <- tibble()

for (t_half in carryover_values) {
  cat(sprintf("Simulating t½ = %.2f weeks... ", t_half))

  iteration_results <- map_dfr(1:n_iterations, function(iter) {
    # Generate data
    data <- generate_olbdc_data(
      n_participants, measurement_weeks, t_half,
      BR_rate, true_beta_bm, noise_sd
    )

    # Analyze WITHOUT carryover modeling
    result_simple <- analyze_data(data, model_carryover = FALSE)
    result_simple$model <- "Without Carryover"

    # Analyze WITH carryover modeling (using 0.2 as assumed halflife)
    result_carryover <- analyze_data(
      data,
      model_carryover = TRUE,
      carryover_halflife = 0.2
    )
    result_carryover$model <- "With Carryover"

    bind_rows(result_simple, result_carryover) |>
      mutate(
        iteration = iter,
        true_t_half = t_half
      )
  })

  results <- bind_rows(results, iteration_results)
  cat("done\n")
}

# ============================================================================
# Summarize Results
# ============================================================================

summary_stats <- results |>
  filter(!is.na(p_value)) |>
  group_by(true_t_half, model) |>
  summarize(
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    power = mean(p_value < 0.05, na.rm = TRUE),
    n_valid = n(),
    .groups = "drop"
  ) |>
  mutate(
    # Approximate expected beta based on true_beta_bm and BR_rate
    expected_beta = -true_beta_bm * BR_rate,  # Negative because outcome decreases
    bias_pct = (mean_beta - expected_beta) / abs(expected_beta) * 100
  )

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESULTS SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Print table
summary_stats |>
  select(true_t_half, model, power, mean_beta, bias_pct) |>
  pivot_wider(
    names_from = model,
    values_from = c(power, mean_beta, bias_pct)
  ) |>
  print(n = 20)

# ============================================================================
# Generate Plots
# ============================================================================

# Power curve plot
p_power <- ggplot(summary_stats, aes(x = true_t_half, y = power, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(breaks = carryover_values) +
  labs(
    title = "Statistical Power Across Carryover Spectrum",
    subtitle = "OL+BDC Design: Non-Monotonic Pattern Without Carryover Modeling",
    x = "True Carryover Half-Life (weeks)",
    y = "Power (at α = 0.05)",
    color = "Analysis Model",
    caption = sprintf(
      "Source: simulation_carryover_spectrum.R | %s | n=%d, iter=%d",
      format(Sys.time(), "%Y-%m-%d"),
      n_participants,
      n_iterations
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  "../output/carryover_spectrum_power.png",
  p_power,
  width = 10,
  height = 6,
  dpi = 150
)

# Bias plot
p_bias <- ggplot(summary_stats, aes(x = true_t_half, y = bias_pct, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(breaks = carryover_values) +
  labs(
    title = "Estimation Bias Across Carryover Spectrum",
    subtitle = "Percent bias in biomarker moderation effect estimate",
    x = "True Carryover Half-Life (weeks)",
    y = "Bias (%)",
    color = "Analysis Model",
    caption = sprintf(
      "Source: simulation_carryover_spectrum.R | %s | n=%d, iter=%d",
      format(Sys.time(), "%Y-%m-%d"),
      n_participants,
      n_iterations
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  "../output/carryover_spectrum_bias.png",
  p_bias,
  width = 10,
  height = 6,
  dpi = 150
)

# Combined plot showing the tradeoff
summary_wide <- summary_stats |>
  filter(model == "Without Carryover") |>
  select(true_t_half, power, bias_pct)

p_tradeoff <- ggplot(summary_wide, aes(x = true_t_half)) +
  geom_line(aes(y = power * 100), color = "#2E86AB", linewidth = 1.2) +
  geom_point(aes(y = power * 100), color = "#2E86AB", size = 3) +
  geom_line(aes(y = -bias_pct), color = "#A23B72", linewidth = 1.2, linetype = "dashed") +
  geom_point(aes(y = -bias_pct), color = "#A23B72", size = 3) +
  scale_x_continuous(breaks = carryover_values) +
  scale_y_continuous(
    name = "Power (%)",
    sec.axis = sec_axis(~ -., name = "Bias (%)")
  ) +
  labs(
    title = "Bias-Variance Tradeoff in Model Misspecification",
    subtitle = "Without Carryover Modeling: Power (blue) vs. Bias (red, dashed)",
    x = "True Carryover Half-Life (weeks)",
    caption = sprintf(
      "Source: simulation_carryover_spectrum.R | %s | n=%d, iter=%d\nNote: Power drops sharply at short carryover, then partially recovers",
      format(Sys.time(), "%Y-%m-%d"),
      n_participants,
      n_iterations
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left = element_text(color = "#2E86AB"),
    axis.title.y.right = element_text(color = "#A23B72"),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  "../output/carryover_spectrum_tradeoff.png",
  p_tradeoff,
  width = 10,
  height = 6,
  dpi = 150
)

# Save results
save(
  results,
  summary_stats,
  file = "../output/carryover_spectrum_results.RData"
)

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Output files saved:\n")
cat("  ../output/carryover_spectrum_power.png\n")
cat("  ../output/carryover_spectrum_bias.png\n")
cat("  ../output/carryover_spectrum_tradeoff.png\n")
cat("  ../output/carryover_spectrum_results.RData\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
