# Debug script: Check if effect size is being recovered correctly
# If model is correctly specified, we should get unbiased estimates
# even if power drops due to reduced precision

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS)
  library(zoo)
})

set.seed(42)

n_participants <- 70
n_sim <- 200
BR_rate <- 0.5
ER_rate <- 0.2
TR_rate <- 0.1
biomarker_mean <- 5.0
biomarker_sd <- 2.0
baseline_mean <- 10.0
bm_mod <- 0.65  # True biomarker moderation effect

# True effect: BR_rate * bm_mod = 0.5 * 0.65 = 0.325
true_effect <- BR_rate * bm_mod
cat("True effect (BR_rate * bm_mod):", true_effect, "\n\n")

measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

# Simplified participant generation (just random effects)
generate_simple_data <- function(n_participants, weeks, bm_mod,
                                  carryover_t_half, model_carryover) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))

  # Create hybrid design
  trial_data <- expand_grid(
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

  # Add participant-level random effects
  participant_info <- tibble(
    participant_id = 1:n_participants,
    biomarker = rnorm(n_participants, biomarker_mean, biomarker_sd),
    baseline = rnorm(n_participants, baseline_mean, 2.0),
    random_intercept = rnorm(n_participants, 0, 2.0)
  )

  trial_data <- trial_data |>
    left_join(participant_info, by = "participant_id") |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      weeks_with_expectancy = cumsum(expectancy),
      weeks_in_trial = week - min(week),
      just_discontinued = treatment == 0 & lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(just_discontinued, week, NA_real_),
      discontinuation_week = zoo::na.locf(discontinuation_week, na.rm = FALSE),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        week - discontinuation_week,
        0
      ),
      bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
      effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered),
      br_on_drug = weeks_on_drug * effective_BR_rate,
      br_at_discont = if_else(just_discontinued, lag(br_on_drug, default = 0), NA_real_),
      br_at_discont = zoo::na.locf(br_at_discont, na.rm = FALSE),
      carryover_effect = if_else(
        treatment == 0 & carryover_t_half > 0,
        (1/2)^(tsd / carryover_t_half),
        0
      ),
      br_off_drug = if_else(
        treatment == 0 & !is.na(br_at_discont) & carryover_t_half > 0,
        br_at_discont * carryover_effect,
        0
      ),
      BR_mean = if_else(treatment == 1, br_on_drug, br_off_drug),
      weeks_at_discont = if_else(just_discontinued, lag(weeks_on_drug, default = 0), NA_real_),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE),
      effective_drug_weeks = if_else(
        treatment == 1,
        weeks_on_drug,
        if_else(!is.na(weeks_at_discont) & carryover_t_half > 0,
                weeks_at_discont * carryover_effect,
                0)
      ),
      ER_mean = weeks_with_expectancy * ER_rate,
      TR_mean = weeks_in_trial * TR_rate,
      noise = rnorm(n(), 0, 1.8),
      response = baseline + random_intercept + BR_mean + ER_mean + TR_mean + noise
    ) |>
    ungroup() |>
    mutate(bm_centered = biomarker - mean(biomarker))

  trial_data
}

# Run simulations for different carryover levels
carryover_levels <- c(0, 0.1, 0.2)  # Hendrickson values (weeks)

results <- list()

for (co in carryover_levels) {
  cat("Testing carryover t_half =", co, "\n")

  estimates <- numeric(n_sim)
  se_values <- numeric(n_sim)
  p_values <- numeric(n_sim)

  for (i in 1:n_sim) {
    set.seed(i * 1000)
    data <- generate_simple_data(n_participants, measurement_weeks,
                                  bm_mod, co, model_carryover = TRUE)

    model <- lmer(
      response ~ effective_drug_weeks * bm_centered + week + (1 | participant_id),
      data = data,
      REML = TRUE
    )

    coefs <- summary(model)$coefficients
    idx <- which(rownames(coefs) == "effective_drug_weeks:bm_centered")

    estimates[i] <- coefs[idx, "Estimate"]
    se_values[i] <- coefs[idx, "Std. Error"]
    t_val <- coefs[idx, "t value"]
    p_values[i] <- 2 * pnorm(-abs(t_val))
  }

  results[[as.character(co)]] <- tibble(
    carryover_t_half = co,
    true_effect = true_effect,
    mean_estimate = mean(estimates),
    sd_estimate = sd(estimates),
    mean_se = mean(se_values),
    bias = mean(estimates) - true_effect,
    bias_pct = 100 * (mean(estimates) - true_effect) / true_effect,
    power = mean(p_values < 0.05),
    coverage_95 = mean(abs(estimates - true_effect) < 1.96 * se_values)
  )

  cat(sprintf("  Mean estimate: %.4f (true: %.4f, bias: %.1f%%)\n",
              mean(estimates), true_effect,
              100 * (mean(estimates) - true_effect) / true_effect))
  cat(sprintf("  Mean SE: %.4f, SD of estimates: %.4f\n",
              mean(se_values), sd(estimates)))
  cat(sprintf("  Power: %.1f%%\n", 100 * mean(p_values < 0.05)))
  cat(sprintf("  95%% CI coverage: %.1f%%\n\n",
              100 * mean(abs(estimates - true_effect) < 1.96 * se_values)))
}

# Combine and display results
results_df <- bind_rows(results)
print(results_df)

cat("\n", strrep("=", 70), "\n")
cat("INTERPRETATION:\n")
cat(strrep("=", 70), "\n")
cat("If bias is near zero across all carryover levels, the model is correctly\n")
cat("specified. Power drops are then due to increased SE (reduced precision),\n")
cat("not model misspecification.\n\n")
cat("If bias increases with carryover, there's a mismatch between DGP and model.\n")
