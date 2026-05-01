#!/usr/bin/env Rscript
# diagnose_nonmonotonic_power.R
#
# Investigate why power is non-monotonic with carryover half-life
# in the OL+BDC design WITHOUT carryover modeling
#
# Hypothesis: With longer carryover, the drug effect persists so long
# that weeks_on_drug (the simple predictor) becomes a better proxy
# for the true cumulative effect

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

set.seed(42)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSING NON-MONOTONIC POWER IN OL+BDC WITHOUT CARRYOVER MODEL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Parameters
n_participants <- 70
BR_rate <- 0.5
biomarker_mean <- 5.0
biomarker_sd <- 2.0

# OL+BDC design: Open-label weeks 1-8, then blinded discontinuation
# Simplified schedule
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

create_olbdc_design <- function(n_participants, weeks) {
  expand_grid(
    participant_id = 1:n_participants,
    week = weeks
  ) |>
    mutate(
      treatment = if_else(week <= 8, 1, 0),
      expectancy = if_else(week <= 8, 1.0, 0.5)
    )
}

# Generate data with different carryover levels
carryover_levels <- c(0, 0.5, 1.0, 2.0)  # Old values from the heatmap

cat("1. PREDICTOR STATISTICS BY CARRYOVER LEVEL\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

for (t_half in carryover_levels) {
  design <- create_olbdc_design(n_participants, measurement_weeks)

  # Add participant info
  participant_info <- tibble(
    participant_id = 1:n_participants,
    biomarker = rnorm(n_participants, biomarker_mean, biomarker_sd)
  )

  data <- design |>
    left_join(participant_info, by = "participant_id") |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
      just_discontinued = treatment == 0 & lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(just_discontinued, week, NA_real_),
      discontinuation_week = zoo::na.locf(discontinuation_week, na.rm = FALSE),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        week - discontinuation_week,
        0
      ),
      weeks_at_discont = if_else(just_discontinued, lag(weeks_on_drug, default = 0), NA_real_),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE),
      # True carryover effect (what's actually happening in DGP)
      carryover_retention = if_else(
        treatment == 0 & t_half > 0,
        (1/2)^(tsd / t_half),
        0
      ),
      # True biological response (what we're trying to detect)
      true_BR = case_when(
        treatment == 1 ~ weeks_on_drug * BR_rate * (1 + 0.5 * bm_centered),
        t_half > 0 ~ weeks_at_discont * BR_rate * (1 + 0.5 * bm_centered) * carryover_retention,
        TRUE ~ 0
      ),
      # Analysis predictor WITHOUT carryover model
      predictor_simple = weeks_on_drug,
      # Analysis predictor WITH carryover model
      predictor_carryover = case_when(
        treatment == 1 ~ weeks_on_drug,
        t_half > 0 ~ weeks_at_discont * carryover_retention,
        TRUE ~ 0
      )
    ) |>
    ungroup()

  # Calculate key statistics
  cat(sprintf("Carryover t½ = %.1f weeks:\n", t_half))

  # Variance of predictors
  var_simple <- var(data$predictor_simple)
  var_carryover <- var(data$predictor_carryover)

  # Correlation between true effect and predictors
  cor_simple_true <- cor(data$predictor_simple * data$bm_centered,
                         data$true_BR, use = "complete.obs")
  cor_carryover_true <- cor(data$predictor_carryover * data$bm_centered,
                            data$true_BR, use = "complete.obs")

  # What fraction of variance in true_BR is explained by simple predictor?
  r2_simple <- cor(data$predictor_simple, data$true_BR, use = "complete.obs")^2
  r2_carryover <- cor(data$predictor_carryover, data$true_BR, use = "complete.obs")^2

  cat(sprintf("  Var(weeks_on_drug): %.2f\n", var_simple))
  cat(sprintf("  Var(effective_drug_weeks): %.2f\n", var_carryover))
  cat(sprintf("  Cor(simple×bm, true_BR): %.3f\n", cor_simple_true))
  cat(sprintf("  Cor(carryover×bm, true_BR): %.3f\n", cor_carryover_true))
  cat(sprintf("  R² simple vs true_BR: %.3f\n", r2_simple))
  cat(sprintf("  R² carryover vs true_BR: %.3f\n", r2_carryover))
  cat("\n")
}

cat("\n2. DETAILED VIEW: OFF-DRUG OBSERVATIONS\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

cat("In OL+BDC, discontinuation happens after week 8.\n")
cat("Off-drug observations are at weeks 9, 10, 11, 12, 16, 20.\n\n")

for (t_half in carryover_levels) {
  cat(sprintf("t½ = %.1f weeks - Carryover retention at off-drug timepoints:\n", t_half))

  off_drug_weeks <- c(9, 10, 11, 12, 16, 20)
  discontinuation_week <- 9  # First off-drug week

  for (w in off_drug_weeks) {
    tsd <- w - discontinuation_week
    if (t_half > 0) {
      retention <- (1/2)^(tsd / t_half)
    } else {
      retention <- 0
    }
    cat(sprintf("  Week %2d (tsd=%2d): retention = %.3f\n", w, tsd, retention))
  }
  cat("\n")
}

cat("\n3. THE KEY INSIGHT\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

cat("When using weeks_on_drug (simple predictor) WITHOUT carryover modeling:\n\n")

cat("At t½=0 (no carryover):\n")
cat("  - Off-drug: true effect = 0, predictor = constant (weeks accumulated)\n")
cat("  - Clean separation between on-drug (effect) and off-drug (no effect)\n")
cat("  - HIGH POWER because contrast is clean\n\n")

cat("At t½=0.5 weeks (short carryover):\n")
cat("  - Off-drug: true effect = decaying, predictor = constant\n")
cat("  - MISMATCH: model assumes no effect, but there IS residual effect\n")
cat("  - This residual effect is UNMODELED NOISE that hurts detection\n")
cat("  - LOWEST POWER due to maximum mismatch\n\n")

cat("At t½=2 weeks (long carryover):\n")
cat("  - Off-drug: true effect decays SLOWLY (still ~50% after 2 weeks)\n")
cat("  - The 'off-drug' periods still have substantial drug effect\n")
cat("  - weeks_on_drug is constant during off-drug, but so is (roughly) the effect!\n")
cat("  - Less mismatch between predictor and true effect\n")
cat("  - POWER RECOVERS somewhat\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CONCLUSION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("The non-monotonicity occurs because:\n\n")
cat("1. At t½=0: Perfect model specification (no carryover in DGP or model)\n\n")
cat("2. At t½=0.5: MAXIMUM misspecification\n")
cat("   - Carryover exists but decays quickly\n")
cat("   - Creates 'noisy' off-drug observations that confuse the model\n")
cat("   - The model sees varying responses during 'off-drug' periods\n")
cat("     but has no predictor to explain them\n\n")
cat("3. At t½=2: Reduced misspecification paradoxically\n")
cat("   - Carryover persists so long that off-drug periods\n")
cat("     look almost like on-drug periods\n")
cat("   - Less contrast, but also less unexplained variance\n")
cat("   - The simple model captures more of the signal\n\n")

cat("This is a form of 'bias-variance tradeoff' in model misspecification.\n")
