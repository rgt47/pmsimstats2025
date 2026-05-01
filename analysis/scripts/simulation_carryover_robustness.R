#!/usr/bin/env Rscript
# simulation_carryover_robustness.R
#
# Compare robustness of carryover-adjusted analysis across different
# data-generating processes (DGPs):
#   - Exponential decay: effect × (1/2)^(tsd / halflife)
#   - Linear decay: effect × max(0, 1 - tsd / (2 × halflife))
#   - Weibull decay: effect × exp(-(tsd / λ)^k)
#
# For each DGP, analyze with and without carryover term to assess:
#   1. Power to detect biomarker moderation
#   2. Bias in effect estimates
#   3. Coverage of 95% confidence intervals
#
# Key question: How robust is exponential carryover modeling when
# the true decay follows a different functional form?

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS)
  library(zoo)
  library(conflicted)
})

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

set.seed(2025)

# ============================================================================
# PARAMETERS
# ============================================================================

n_participants <- 70
n_iterations <- 200

# Three-factor response model - RATE-BASED (points per week)
BR_rate <- 0.5   # Biological Response: drug improvement rate
ER_rate <- 0.2   # Expectancy Response: placebo improvement rate
TR_rate <- 0.1   # Time-variant Response: natural improvement rate

baseline_mean <- 10.0
between_subject_sd <- 2.0
within_subject_sd <- 1.8
biomarker_mean <- 5.0
biomarker_sd <- 2.0

# Measurement schedule - hybrid design clustered timepoints
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CARRYOVER DECAY MODEL ROBUSTNESS ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Simulation parameters:\n")
cat(sprintf("  Participants: %d\n", n_participants))
cat(sprintf("  Iterations: %d\n", n_iterations))
cat(sprintf("  BR rate: %.2f, ER rate: %.2f, TR rate: %.2f\n",
            BR_rate, ER_rate, TR_rate))
cat(sprintf("  Measurement weeks: %s\n",
            paste(measurement_weeks, collapse = ", ")))
cat("\n")

# ============================================================================
# CARRYOVER DECAY FUNCTIONS
# ============================================================================

#' Exponential decay carryover
#' Standard model: effect decays with half-life
#' @param tsd Time since discontinuation (weeks)
#' @param halflife Half-life parameter (weeks)
#' @param effect_at_discont Effect magnitude at discontinuation
carryover_exponential <- function(tsd, halflife, effect_at_discont) {

  if (halflife <= 0 || tsd <= 0) return(0)
  effect_at_discont * (1/2)^(tsd / halflife)
}

#' Linear decay carryover
#' Effect decays linearly to zero over 2 × halflife
#' Matched so that at t = halflife, retention = 0.5 (same as exponential)
#' @param tsd Time since discontinuation (weeks)
#' @param halflife Time to reach 50% (weeks), full decay at 2 × halflife
#' @param effect_at_discont Effect magnitude at discontinuation
carryover_linear <- function(tsd, halflife, effect_at_discont) {
  if (halflife <= 0 || tsd <= 0) return(0)
  total_decay_time <- 2 * halflife
  retention <- pmax(0, 1 - tsd / total_decay_time)
  effect_at_discont * retention
}

#' Weibull decay carryover
#' Flexible decay: k < 1 = slow then fast, k > 1 = fast then slow
#' Parameterized so halflife matches exponential at t = halflife
#' @param tsd Time since discontinuation (weeks)
#' @param halflife Effective half-life (weeks)
#' @param effect_at_discont Effect magnitude at discontinuation
#' @param k Shape parameter (k = 1 gives exponential)
carryover_weibull <- function(tsd, halflife, effect_at_discont, k = 1.5) {
  if (halflife <= 0 || tsd <= 0) return(0)
  lambda <- halflife / (log(2)^(1/k))
  effect_at_discont * exp(-(tsd / lambda)^k)
}

#' Master carryover function
#' @param tsd Time since discontinuation
#' @param halflife Half-life parameter
#' @param effect_at_discont Effect at discontinuation
#' @param decay_type One of "exponential", "linear", "weibull"
#' @param weibull_k Shape parameter for Weibull (default 1.5)
calculate_carryover <- function(tsd, halflife, effect_at_discont,
                                 decay_type = "exponential",
                                 weibull_k = 1.5) {
  if (is.na(tsd) || tsd <= 0 || halflife <= 0) return(0)

  switch(decay_type,
    "exponential" = carryover_exponential(tsd, halflife, effect_at_discont),
    "linear" = carryover_linear(tsd, halflife, effect_at_discont),
    "weibull" = carryover_weibull(tsd, halflife, effect_at_discont, weibull_k),
    stop("Unknown decay type: ", decay_type)
  )
}

# Vectorized version
calculate_carryover_vec <- Vectorize(calculate_carryover)

# ============================================================================
# VISUALIZE DECAY CURVES
# ============================================================================

cat("Carryover decay models (halflife = 0.2 weeks):\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")

tsd_seq <- seq(0, 1, by = 0.05)
halflife_demo <- 0.2

decay_comparison <- tibble(
  tsd = rep(tsd_seq, 3),
  decay_type = rep(c("exponential", "linear", "weibull"), each = length(tsd_seq)),
  retention = c(
    sapply(tsd_seq, \(t) if(t == 0) 1 else carryover_exponential(t, halflife_demo, 1)),
    sapply(tsd_seq, \(t) if(t == 0) 1 else carryover_linear(t, halflife_demo, 1)),
    sapply(tsd_seq, \(t) if(t == 0) 1 else carryover_weibull(t, halflife_demo, 1, k = 1.5))
  )
)

cat("\nRetention at key timepoints (halflife = 0.2 weeks):\n")
decay_comparison |>
  filter(tsd %in% c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) |>
  pivot_wider(names_from = decay_type, values_from = retention) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print()

cat("\nNote: All models calibrated so retention ≈ 0.5 at t = halflife\n")
cat("  - Exponential: smooth geometric decay\
")
cat("  - Linear: constant rate to zero at 2×halflife\n")
cat("  - Weibull (k=1.5): initially faster, then slower decay\n\n")

# ============================================================================
# PARAMETER GRID
# ============================================================================

param_grid <- expand_grid(
  design = c("hybrid", "crossover"),
  decay_type = c("exponential", "linear", "weibull"),
  biomarker_moderation = c(0.3, 0.6),
  carryover_halflife = c(0, 0.1, 0.2),  # Hendrickson values (weeks)
  model_carryover = c(TRUE, FALSE)
) |>
  filter(!(carryover_halflife == 0 & model_carryover == FALSE)) |>
  mutate(condition_id = row_number())

cat("Parameter grid:\n")
cat(sprintf("  Designs: %s\n", paste(unique(param_grid$design), collapse = ", ")))
cat(sprintf("  Decay types: %s\n", paste(unique(param_grid$decay_type), collapse = ", ")))
cat(sprintf("  Biomarker moderation: %s\n",
            paste(unique(param_grid$biomarker_moderation), collapse = ", ")))
cat(sprintf("  Carryover half-lives: %s weeks\n",
            paste(unique(param_grid$carryover_halflife), collapse = ", ")))
cat(sprintf("  Model carryover: TRUE and FALSE\n"))
cat(sprintf("  Total conditions: %d\n\n", nrow(param_grid)))

# ============================================================================
# TRIAL DESIGN TEMPLATES
# ============================================================================

create_hybrid_design <- function(n_participants, weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = weeks
  ) |>
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        week %in% c(4, 8) ~ 1,
        week == 9 ~ 1,
        week == 10 & path %in% c(1, 2) ~ 1,
        week == 10 & path %in% c(3, 4) ~ 0,
        week %in% c(11, 12) ~ 0,
        week == 16 & path %in% c(1, 3) ~ 1,
        week == 16 & path %in% c(2, 4) ~ 0,
        week == 20 & path %in% c(1, 3) ~ 0,
        week == 20 & path %in% c(2, 4) ~ 1,
        TRUE ~ NA_real_
      ),
      expectancy = if_else(week %in% c(4, 8), 1.0, 0.5)
    )
}

create_crossover_design <- function(n_participants, weeks) {
  sequence_assignment <- sample(rep(1:2, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = weeks
  ) |>
    mutate(
      sequence = sequence_assignment[participant_id],
      treatment = case_when(
        week %in% c(4, 8) & sequence == 1 ~ 1,
        week %in% c(4, 8) & sequence == 2 ~ 0,
        week %in% c(9, 10, 11, 12) & sequence == 1 ~ 0,
        week %in% c(9, 10, 11, 12) & sequence == 2 ~ 1,
        week %in% c(16, 20) & sequence == 1 ~ 1,
        week %in% c(16, 20) & sequence == 2 ~ 0,
        TRUE ~ NA_real_
      ),
      expectancy = 0.5
    )
}

# ============================================================================
# DATA GENERATION FUNCTION
# ============================================================================

generate_trial_data <- function(n_participants, weeks, bm_mod,
                                 carryover_halflife, decay_type) {
  design_data <- create_hybrid_design(n_participants, weeks)

  participant_info <- tibble(
    participant_id = 1:n_participants,
    biomarker = rnorm(n_participants, biomarker_mean, biomarker_sd),
    baseline = rnorm(n_participants, baseline_mean, between_subject_sd),
    random_intercept = rnorm(n_participants, 0, between_subject_sd)
  )

  trial_data <- design_data |>
    left_join(participant_info, by = "participant_id") |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      weeks_with_expectancy = cumsum(expectancy),
      weeks_in_trial = week - min(week),
      bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
      effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered),
      br_on_drug = weeks_on_drug * effective_BR_rate,
      just_discontinued = treatment == 0 & lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(just_discontinued, week, NA_real_),
      discontinuation_week = zoo::na.locf(discontinuation_week, na.rm = FALSE),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        week - discontinuation_week,
        0
      ),
      br_at_discont = if_else(
        just_discontinued, lag(br_on_drug, default = 0), NA_real_
      ),
      br_at_discont = zoo::na.locf(br_at_discont, na.rm = FALSE)
    ) |>
    ungroup()

  if (carryover_halflife > 0) {
    trial_data <- trial_data |>
      rowwise() |>
      mutate(
        carryover_effect = calculate_carryover(
          tsd, carryover_halflife, br_at_discont, decay_type
        )
      ) |>
      ungroup()
  } else {
    trial_data <- trial_data |>
      mutate(carryover_effect = 0)
  }

  trial_data <- trial_data |>
    mutate(
      BR_mean = if_else(treatment == 1, br_on_drug, carryover_effect),
      ER_mean = weeks_with_expectancy * ER_rate,
      TR_mean = weeks_in_trial * TR_rate,
      noise = rnorm(n(), 0, within_subject_sd),
      response = baseline + random_intercept + BR_mean + ER_mean + TR_mean + noise
    ) |>
    group_by(participant_id) |>
    mutate(
      weeks_at_discont = if_else(
        just_discontinued, lag(weeks_on_drug, default = 0), NA_real_
      ),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE),
      effective_drug_weeks = if_else(
        treatment == 1,
        weeks_on_drug,
        if_else(
          !is.na(weeks_at_discont) & carryover_halflife > 0,
          weeks_at_discont * (1/2)^(tsd / carryover_halflife),
          0
        )
      )
    ) |>
    ungroup()

  # Use SAME standardization as in data generation for analysis
  # bm_centered is already standardized during data generation
  # Keep it consistent for analysis

  trial_data

  trial_data
}

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

fit_model_with_carryover <- function(data) {
  tryCatch({
    model <- lmer(
      response ~ effective_drug_weeks * bm_centered + week + (1 | participant_id),
      data = data,
      REML = TRUE,
      control = lmerControl(optimizer = "bobyqa")
    )

    coefs <- summary(model)$coefficients
    idx <- which(rownames(coefs) == "effective_drug_weeks:bm_centered")

    if (length(idx) == 0) {
      return(list(estimate = NA, se = NA, p_value = NA, converged = FALSE))
    }

    estimate <- coefs[idx, "Estimate"]
    se <- coefs[idx, "Std. Error"]
    t_val <- coefs[idx, "t value"]
    p_value <- 2 * pnorm(-abs(t_val))

    list(estimate = estimate, se = se, p_value = p_value, converged = TRUE)
  }, error = function(e) {
    list(estimate = NA, se = NA, p_value = NA, converged = FALSE)
  })
}

fit_model_without_carryover <- function(data) {
  tryCatch({
    model <- lmer(
      response ~ weeks_on_drug * bm_centered + week + (1 | participant_id),
      data = data,
      REML = TRUE,
      control = lmerControl(optimizer = "bobyqa")
    )

    coefs <- summary(model)$coefficients
    idx <- which(rownames(coefs) == "weeks_on_drug:bm_centered")

    if (length(idx) == 0) {
      return(list(estimate = NA, se = NA, p_value = NA, converged = FALSE))
    }

    estimate <- coefs[idx, "Estimate"]
    se <- coefs[idx, "Std. Error"]
    t_val <- coefs[idx, "t value"]
    p_value <- 2 * pnorm(-abs(t_val))

    list(estimate = estimate, se = se, p_value = p_value, converged = TRUE)
  }, error = function(e) {
    list(estimate = NA, se = NA, p_value = NA, converged = FALSE)
  })
}

# ============================================================================
# RUN SIMULATION
# ============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RUNNING SIMULATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

results_list <- list()

for (i in seq_len(nrow(param_grid))) {
  params <- param_grid[i, ]

  cat(sprintf(
    "\rCondition %d/%d: %s design, %s decay, bm=%.1f, t½=%.1f, model=%s",
    i, nrow(param_grid),
    params$design, params$decay_type,
    params$biomarker_moderation, params$carryover_halflife,
    params$model_carryover
  ))

  true_effect <- BR_rate * params$biomarker_moderation

  estimates <- numeric(n_iterations)
  se_values <- numeric(n_iterations)
  p_values <- numeric(n_iterations)
  converged <- logical(n_iterations)

  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    data <- generate_trial_data(
      n_participants,
      measurement_weeks,
      params$biomarker_moderation,
      params$carryover_halflife,
      params$decay_type
    )

    if (params$model_carryover) {
      result <- fit_model_with_carryover(data)
    } else {
      result <- fit_model_without_carryover(data)
    }

    estimates[iter] <- result$estimate
    se_values[iter] <- result$se
    p_values[iter] <- result$p_value
    converged[iter] <- result$converged
  }

  valid <- converged & !is.na(estimates)

  results_list[[i]] <- tibble(
    condition_id = params$condition_id,
    design = params$design,
    decay_type = params$decay_type,
    biomarker_moderation = params$biomarker_moderation,
    carryover_halflife = params$carryover_halflife,
    model_carryover = params$model_carryover,
    true_effect = true_effect,
    mean_estimate = mean(estimates[valid], na.rm = TRUE),
    sd_estimate = sd(estimates[valid], na.rm = TRUE),
    mean_se = mean(se_values[valid], na.rm = TRUE),
    bias = mean(estimates[valid], na.rm = TRUE) - true_effect,
    bias_pct = 100 * (mean(estimates[valid], na.rm = TRUE) - true_effect) /
               true_effect,
    rmse = sqrt(mean((estimates[valid] - true_effect)^2, na.rm = TRUE)),
    power = mean(p_values[valid] < 0.05, na.rm = TRUE),
    coverage_95 = mean(
      abs(estimates[valid] - true_effect) < 1.96 * se_values[valid],
      na.rm = TRUE
    ),
    convergence_rate = mean(converged),
    n_valid = sum(valid)
  )
}

cat("\n\n")

results <- bind_rows(results_list)

# ============================================================================
# SUMMARY RESULTS
# ============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESULTS SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("1. POWER BY DECAY TYPE AND ANALYSIS MODEL\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

power_summary <- results |>
  filter(carryover_halflife > 0) |>
  group_by(decay_type, model_carryover, carryover_halflife) |>
  summarise(
    mean_power = mean(power),
    mean_bias_pct = mean(bias_pct),
    mean_coverage = mean(coverage_95),
    .groups = "drop"
  ) |>
  arrange(decay_type, carryover_halflife, desc(model_carryover))

print(power_summary, n = 30)

cat("\n2. ROBUSTNESS: BIAS WHEN ANALYSIS MODEL MISSPECIFIED\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("(Data generated with one decay type, analyzed assuming exponential)\n\n")

robustness_summary <- results |>
  filter(carryover_halflife > 0, model_carryover == TRUE) |>
  select(decay_type, carryover_halflife, biomarker_moderation,
         bias_pct, power, coverage_95) |>
  arrange(decay_type, carryover_halflife, biomarker_moderation)

print(robustness_summary, n = 30)

cat("\n3. COMPARISON: WITH vs WITHOUT CARRYOVER MODELING\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

comparison <- results |>
  filter(carryover_halflife > 0) |>
  group_by(decay_type, carryover_halflife, biomarker_moderation, model_carryover) |>
  summarise(
    power = mean(power),
    bias_pct = mean(bias_pct),
    coverage_95 = mean(coverage_95),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = model_carryover,
    values_from = c(power, bias_pct, coverage_95),
    names_glue = "{.value}_{ifelse(model_carryover, 'with', 'without')}"
  ) |>
  mutate(
    power_gain = power_with - power_without,
    bias_reduction = abs(bias_pct_without) - abs(bias_pct_with)
  )

print(comparison |> select(decay_type, carryover_halflife, biomarker_moderation,
                           power_with, power_without, power_gain), n = 30)

cat("\n4. DETAILED RESULTS BY CONDITION\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

results |>
  select(decay_type, carryover_halflife, biomarker_moderation,
         model_carryover, power, bias_pct, coverage_95, rmse) |>
  arrange(decay_type, carryover_halflife, biomarker_moderation,
          desc(model_carryover)) |>
  print(n = 50)

# ============================================================================
# KEY FINDINGS
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("KEY FINDINGS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

avg_power_with <- results |>
  filter(carryover_halflife > 0, model_carryover == TRUE) |>
  pull(power) |>
  mean()

avg_power_without <- results |>
  filter(carryover_halflife > 0, model_carryover == FALSE) |>
  pull(power) |>
  mean()

avg_bias_with <- results |>
  filter(carryover_halflife > 0, model_carryover == TRUE) |>
  pull(bias_pct) |>
  abs() |>
  mean()

avg_bias_without <- results |>
  filter(carryover_halflife > 0, model_carryover == FALSE) |>
  pull(bias_pct) |>
  abs() |>
  mean()

cat(sprintf("Average power WITH carryover modeling: %.1f%%\n",
            100 * avg_power_with))
cat(sprintf("Average power WITHOUT carryover modeling: %.1f%%\n",
            100 * avg_power_without))
cat(sprintf("Power improvement: +%.1f percentage points\n\n",
            100 * (avg_power_with - avg_power_without)))

cat(sprintf("Average |bias| WITH carryover modeling: %.1f%%\n", avg_bias_with))
cat(sprintf("Average |bias| WITHOUT carryover modeling: %.1f%%\n",
            avg_bias_without))

cat("\nROBUSTNESS ACROSS DECAY TYPES:\n")
for (dt in c("exponential", "linear", "weibull")) {
  dt_results <- results |>
    filter(decay_type == dt, carryover_halflife > 0, model_carryover == TRUE)
  cat(sprintf("  %s: mean power = %.1f%%, mean |bias| = %.1f%%\n",
              dt,
              100 * mean(dt_results$power),
              mean(abs(dt_results$bias_pct))))
}

cat("\n")
cat("INTERPRETATION:\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("If power and bias are similar across decay types when using the\n")
cat("exponential carryover model in analysis, the method is ROBUST to\n")
cat("misspecification of the true decay functional form.\n\n")
cat("If linear/Weibull DGPs show worse performance with exponential\n")
cat("analysis model, the method is SENSITIVE to decay specification.\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

save(results, param_grid, file = "../output/carryover_robustness_results.RData")
cat("\nResults saved to: ../output/carryover_robustness_results.RData\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GENERATING VISUALIZATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

plot_data <- results |>
  filter(carryover_halflife > 0) |>
  mutate(
    analysis_label = if_else(model_carryover,
                              "With Carryover Model",
                              "Without Carryover Model"),
    decay_label = factor(decay_type,
                         levels = c("exponential", "linear", "weibull"),
                         labels = c("Exponential", "Linear", "Weibull")),
    halflife_label = sprintf("t½ = %.1f wk", carryover_halflife)
  )

p1 <- ggplot(plot_data,
             aes(x = decay_label, y = power * 100, fill = analysis_label)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_grid(biomarker_moderation ~ halflife_label,
             labeller = labeller(
               biomarker_moderation = \(x) sprintf("BM mod = %.1f", as.numeric(x))
             )) +
  scale_fill_manual(values = c("With Carryover Model" = "#2166AC",
                                "Without Carryover Model" = "#B2182B")) +
  labs(
    title = "Statistical Power by Data-Generating Decay Type",
    subtitle = "Robustness of exponential carryover analysis across different true decay forms",
    x = "True Decay Type (Data Generation)",
    y = "Power (%)",
    fill = "Analysis Model",
    caption = sprintf("Source: simulation_carryover_robustness.R | %s | n=%d, iter=%d",
                      format(Sys.time(), "%Y-%m-%d"), n_participants, n_iterations)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 7, hjust = 1, color = "gray50")
  ) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray50", linewidth = 0.5)

p2 <- ggplot(plot_data,
             aes(x = decay_label, y = bias_pct, fill = analysis_label)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_grid(biomarker_moderation ~ halflife_label,
             labeller = labeller(
               biomarker_moderation = \(x) sprintf("BM mod = %.1f", as.numeric(x))
             )) +
  scale_fill_manual(values = c("With Carryover Model" = "#2166AC",
                                "Without Carryover Model" = "#B2182B")) +
  labs(
    title = "Estimation Bias by Data-Generating Decay Type",
    subtitle = "Percent bias in biomarker moderation effect estimate",
    x = "True Decay Type (Data Generation)",
    y = "Bias (%)",
    fill = "Analysis Model",
    caption = sprintf("Source: simulation_carryover_robustness.R | %s | n=%d, iter=%d",
                      format(Sys.time(), "%Y-%m-%d"), n_participants, n_iterations)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 7, hjust = 1, color = "gray50")
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray30", linewidth = 0.5)

p3 <- ggplot(plot_data,
             aes(x = decay_label, y = coverage_95 * 100, fill = analysis_label)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_grid(biomarker_moderation ~ halflife_label,
             labeller = labeller(
               biomarker_moderation = \(x) sprintf("BM mod = %.1f", as.numeric(x))
             )) +
  scale_fill_manual(values = c("With Carryover Model" = "#2166AC",
                                "Without Carryover Model" = "#B2182B")) +
  labs(
    title = "95% CI Coverage by Data-Generating Decay Type",
    subtitle = "Proportion of simulations where true effect falls within estimated 95% CI",
    x = "True Decay Type (Data Generation)",
    y = "Coverage (%)",
    fill = "Analysis Model",
    caption = sprintf("Source: simulation_carryover_robustness.R | %s | n=%d, iter=%d",
                      format(Sys.time(), "%Y-%m-%d"), n_participants, n_iterations)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 7, hjust = 1, color = "gray50")
  ) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray50", linewidth = 0.5)

ggsave("../output/carryover_robustness_power.png", p1,
       width = 10, height = 7, dpi = 150)
ggsave("../output/carryover_robustness_bias.png", p2,
       width = 10, height = 7, dpi = 150)
ggsave("../output/carryover_robustness_coverage.png", p3,
       width = 10, height = 7, dpi = 150)

cat("Plots saved to:\n")
cat("  ../output/carryover_robustness_power.png\n")
cat("  ../output/carryover_robustness_bias.png\n")
cat("  ../output/carryover_robustness_coverage.png\n")

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SIMULATION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
