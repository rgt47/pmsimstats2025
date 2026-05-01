# Model Specification Diagnostic
# Check whether the analysis model correctly matches the data generating process

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS)
  library(zoo)
})

set.seed(42)

# Parameters matching simulation
n_participants <- 70
BR_rate <- 0.5
ER_rate <- 0.2
TR_rate <- 0.1
biomarker_mean <- 5.0
biomarker_sd <- 2.0
baseline_mean <- 10.0
bm_mod <- 0.65

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL SPECIFICATION DIAGNOSTIC\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Data Generating Process (DGP):\n")
cat("  BR_mean = effective_BR_rate * weeks_on_drug  (on drug)\n")
cat("  BR_mean = br_at_discont * carryover_effect   (off drug)\n")
cat("  where effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered)\n")
cat("  and   br_at_discont = weeks_at_discont * effective_BR_rate\n\n")

cat("Expanding the DGP:\n")
cat("  On drug:  BR_mean = weeks_on_drug * BR_rate * (1 + bm_mod * bm_centered)\n")
cat("  Off drug: BR_mean = weeks_at_discont * BR_rate * (1 + bm_mod * bm_centered) *",
    "carryover_effect\n\n")

cat("This simplifies to:\n")
cat("  BR_mean = effective_drug_weeks * BR_rate * (1 + bm_mod * bm_centered)\n")
cat("  where effective_drug_weeks = weeks_on_drug (on) or",
    "weeks_at_discont * carryover_effect (off)\n\n")

cat("Expanding further:\n")
cat("  BR_mean = effective_drug_weeks * BR_rate +",
    "effective_drug_weeks * BR_rate * bm_mod * bm_centered\n\n")

cat("Therefore, the TRUE coefficients in a linear model are:\n")
cat("  - Coefficient on effective_drug_weeks: BR_rate =", BR_rate, "\n")
cat("  - Coefficient on effective_drug_weeks:bm_centered: BR_rate * bm_mod =",
    BR_rate * bm_mod, "\n\n")

# Generate one dataset to verify
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)
carryover_t_half <- 1.0

path_assignment <- sample(rep(1:4, length.out = n_participants))

trial_data <- expand_grid(
  participant_id = 1:n_participants,
  week = measurement_weeks
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

cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("VERIFICATION: Does effective_drug_weeks match the DGP?\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# Check: BR_mean should equal effective_drug_weeks * effective_BR_rate
trial_data <- trial_data |>
  mutate(
    BR_reconstructed = effective_drug_weeks * effective_BR_rate,
    match = abs(BR_mean - BR_reconstructed) < 1e-10
  )

cat("All BR_mean values match effective_drug_weeks * effective_BR_rate:",
    all(trial_data$match), "\n\n")

# Show a few examples
cat("Sample data (first participant, paths with discontinuation):\n")
trial_data |>
  filter(participant_id == 1) |>
  dplyr::select(week, treatment, weeks_on_drug, tsd, carryover_effect,
         effective_drug_weeks, BR_mean, BR_reconstructed, match) |>
  print(n = 8)

cat("\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL FIT COMPARISON\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# Fit the current model
model_current <- lmer(
  response ~ effective_drug_weeks * bm_centered + week + (1 | participant_id),
  data = trial_data,
  REML = TRUE
)

cat("Current Model: response ~ effective_drug_weeks * bm_centered + week\n")
cat("True coefficient on effective_drug_weeks:bm_centered:", BR_rate * bm_mod, "\n\n")
print(summary(model_current)$coefficients)

cat("\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("ALTERNATIVE: SEPARATE TREATMENT AND CARRYOVER TERMS\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# Create separate treatment and carryover terms
trial_data <- trial_data |>
  mutate(
    carryover_weeks = if_else(
      treatment == 0 & !is.na(weeks_at_discont) & carryover_t_half > 0,
      weeks_at_discont * carryover_effect,
      0
    )
  )

# Fit alternative model with separate terms
model_separate <- lmer(
  response ~ weeks_on_drug * bm_centered + carryover_weeks * bm_centered +
    week + (1 | participant_id),
  data = trial_data,
  REML = TRUE
)

cat("Alternative Model: response ~ weeks_on_drug * bm_centered +",
    "carryover_weeks * bm_centered + week\n")
cat("True coefficient on weeks_on_drug:bm_centered:", BR_rate * bm_mod, "\n")
cat("True coefficient on carryover_weeks:bm_centered:", BR_rate * bm_mod,
    "(same as active treatment)\n\n")
print(summary(model_separate)$coefficients)

cat("\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("KEY INSIGHT\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("The DGP creates BR_mean as:\n")
cat("  BR_mean = effective_drug_weeks * BR_rate * (1 + bm_mod * bm_centered)\n\n")

cat("This is a MULTIPLICATIVE model where biomarker moderation applies to the\n")
cat("RATE of response, not as a separate additive interaction.\n\n")

cat("The linear model approximates this as:\n")
cat("  E[response] = b0 + b1*effective_drug_weeks + b2*bm_centered +\n")
cat("                b3*effective_drug_weeks:bm_centered + ...\n\n")

cat("For this to work correctly:\n")
cat("  - b1 should estimate BR_rate =", BR_rate, "\n")
cat("  - b3 should estimate BR_rate * bm_mod =", BR_rate * bm_mod, "\n\n")

# Check actual estimates
coef_edw <- fixef(model_current)["effective_drug_weeks"]
coef_int <- fixef(model_current)["effective_drug_weeks:bm_centered"]

cat("Actual estimates from current model:\n")
cat("  b1 (effective_drug_weeks):", round(coef_edw, 4), "\n")
cat("  b3 (effective_drug_weeks:bm_centered):", round(coef_int, 4), "\n\n")

cat("Bias:\n")
cat("  b1 bias:", round(coef_edw - BR_rate, 4),
    sprintf("(%.1f%%)", 100 * (coef_edw - BR_rate) / BR_rate), "\n")
cat("  b3 bias:", round(coef_int - BR_rate * bm_mod, 4),
    sprintf("(%.1f%%)", 100 * (coef_int - BR_rate * bm_mod) / (BR_rate * bm_mod)), "\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ISSUE IDENTIFIED: bm_centered is re-computed using sample mean\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("In the DGP, bm_centered = (biomarker - biomarker_mean) / biomarker_sd\n")
cat("But after data generation, bm_centered is overwritten as:\n")
cat("  bm_centered = biomarker - mean(biomarker)  # sample mean, not standardized!\n\n")

cat("This creates a mismatch. The DGP uses standardized biomarker (SD=1),\n")
cat("but the analysis uses unstandardized (SD ~ biomarker_sd).\n\n")

# Check the mismatch
cat("Biomarker summary:\n")
cat("  True mean:", biomarker_mean, ", True SD:", biomarker_sd, "\n")
cat("  Sample mean:", round(mean(trial_data$biomarker), 2),
    ", Sample SD:", round(sd(trial_data$biomarker), 2), "\n")
cat("  bm_centered SD (after re-centering):", round(sd(trial_data$bm_centered), 2), "\n\n")

cat("Expected coefficient on interaction term:\n")
cat("  If using standardized bm_centered: BR_rate * bm_mod =", BR_rate * bm_mod, "\n")
cat("  If using unstandardized bm_centered: BR_rate * bm_mod / biomarker_sd =",
    BR_rate * bm_mod / biomarker_sd, "\n\n")

cat("This scaling difference explains why the estimate might differ from the\n")
cat("expected value based on the true parameters.\n")
