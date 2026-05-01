# Model Specification Diagnostic - Part 2
# Fix the bm_centered scaling issue and verify model recovery

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS)
  library(zoo)
})

set.seed(42)

n_participants <- 70
BR_rate <- 0.5
ER_rate <- 0.2
TR_rate <- 0.1
biomarker_mean <- 5.0
biomarker_sd <- 2.0
baseline_mean <- 10.0
bm_mod <- 0.65

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL SPECIFICATION DIAGNOSTIC - PART 2\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

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

# KEY: Store the standardized bm_centered used in DGP
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
    # KEEP the standardized version used in DGP
    bm_centered_dgp = (biomarker - biomarker_mean) / biomarker_sd,
    effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered_dgp),
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
  ungroup()

# Create BOTH versions of bm_centered for comparison
trial_data <- trial_data |>
  mutate(
    # Version 1: Sample-centered (NOT standardized) - what simulation currently does
    bm_centered_sample = biomarker - mean(biomarker),
    # Version 2: Standardized using sample stats
    bm_centered_std = (biomarker - mean(biomarker)) / sd(biomarker)
  )

cat("TESTING THREE DIFFERENT bm_centered SPECIFICATIONS\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# Model 1: Using sample-centered (current simulation approach)
model_sample <- lmer(
  response ~ effective_drug_weeks * bm_centered_sample + week + (1 | participant_id),
  data = trial_data,
  REML = TRUE
)

# Model 2: Using sample-standardized
model_std <- lmer(
  response ~ effective_drug_weeks * bm_centered_std + week + (1 | participant_id),
  data = trial_data,
  REML = TRUE
)

# Model 3: Using DGP-standardized (should give best estimates)
model_dgp <- lmer(
  response ~ effective_drug_weeks * bm_centered_dgp + week + (1 | participant_id),
  data = trial_data,
  REML = TRUE
)

cat("True values:\n")
cat("  BR_rate:", BR_rate, "\n")
cat("  BR_rate * bm_mod:", BR_rate * bm_mod, "\n\n")

cat("Model 1: bm_centered = biomarker - mean(biomarker) [CURRENT]\n")
cat("  Expected interaction coef: BR_rate * bm_mod / biomarker_sd =",
    BR_rate * bm_mod / biomarker_sd, "\n")
coef1 <- fixef(model_sample)
cat("  Estimated interaction coef:", round(coef1["effective_drug_weeks:bm_centered_sample"], 4), "\n\n")

cat("Model 2: bm_centered = (biomarker - mean) / sd [STANDARDIZED]\n")
cat("  Expected interaction coef: BR_rate * bm_mod =", BR_rate * bm_mod, "\n")
coef2 <- fixef(model_std)
cat("  Estimated interaction coef:", round(coef2["effective_drug_weeks:bm_centered_std"], 4), "\n\n")

cat("Model 3: bm_centered = (biomarker - true_mean) / true_sd [DGP VERSION]\n")
cat("  Expected interaction coef: BR_rate * bm_mod =", BR_rate * bm_mod, "\n")
coef3 <- fixef(model_dgp)
cat("  Estimated interaction coef:", round(coef3["effective_drug_weeks:bm_centered_dgp"], 4), "\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SIMULATION WITH CORRECT SPECIFICATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Run multiple iterations to check unbiasedness
n_sim <- 100
estimates <- matrix(NA, nrow = n_sim, ncol = 4)
colnames(estimates) <- c("edw_coef", "interaction_coef", "edw_pval", "interaction_pval")

cat("Running", n_sim, "simulations with standardized bm_centered...\n")

for (i in 1:n_sim) {
  set.seed(i * 1000)

  path_assignment <- sample(rep(1:4, length.out = n_participants))

  sim_data <- expand_grid(
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

  sim_participant_info <- tibble(
    participant_id = 1:n_participants,
    biomarker = rnorm(n_participants, biomarker_mean, biomarker_sd),
    baseline = rnorm(n_participants, baseline_mean, 2.0),
    random_intercept = rnorm(n_participants, 0, 2.0)
  )

  sim_data <- sim_data |>
    left_join(sim_participant_info, by = "participant_id") |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      just_discontinued = treatment == 0 & lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(just_discontinued, week, NA_real_),
      discontinuation_week = zoo::na.locf(discontinuation_week, na.rm = FALSE),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        week - discontinuation_week,
        0
      ),
      bm_centered_dgp = (biomarker - biomarker_mean) / biomarker_sd,
      effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered_dgp),
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
      weeks_with_expectancy = cumsum(expectancy),
      weeks_in_trial = week - min(week),
      ER_mean = weeks_with_expectancy * ER_rate,
      TR_mean = weeks_in_trial * TR_rate,
      noise = rnorm(n(), 0, 1.8),
      response = baseline + random_intercept + BR_mean + ER_mean + TR_mean + noise
    ) |>
    ungroup()

  # Use STANDARDIZED bm_centered for analysis (but using sample stats)
  sim_data <- sim_data |>
    mutate(bm_centered = (biomarker - mean(biomarker)) / sd(biomarker))

  model <- lmer(
    response ~ effective_drug_weeks * bm_centered + week + (1 | participant_id),
    data = sim_data,
    REML = TRUE
  )

  coefs <- summary(model)$coefficients
  estimates[i, 1] <- coefs["effective_drug_weeks", "Estimate"]
  estimates[i, 2] <- coefs["effective_drug_weeks:bm_centered", "Estimate"]

  t_edw <- coefs["effective_drug_weeks", "t value"]
  t_int <- coefs["effective_drug_weeks:bm_centered", "t value"]
  estimates[i, 3] <- 2 * pnorm(-abs(t_edw))
  estimates[i, 4] <- 2 * pnorm(-abs(t_int))
}

cat("\nResults with STANDARDIZED bm_centered:\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("effective_drug_weeks coefficient:\n")
cat("  True value:", BR_rate, "\n")
cat("  Mean estimate:", round(mean(estimates[, 1]), 4), "\n")
cat("  SD of estimates:", round(sd(estimates[, 1]), 4), "\n")
cat("  Bias:", round(mean(estimates[, 1]) - BR_rate, 4),
    sprintf("(%.1f%%)", 100 * (mean(estimates[, 1]) - BR_rate) / BR_rate), "\n\n")

cat("effective_drug_weeks:bm_centered coefficient:\n")
cat("  True value:", BR_rate * bm_mod, "\n")
cat("  Mean estimate:", round(mean(estimates[, 2]), 4), "\n")
cat("  SD of estimates:", round(sd(estimates[, 2]), 4), "\n")
cat("  Bias:", round(mean(estimates[, 2]) - BR_rate * bm_mod, 4),
    sprintf("(%.1f%%)", 100 * (mean(estimates[, 2]) - BR_rate * bm_mod) / (BR_rate * bm_mod)), "\n\n")

cat("Power (p < 0.05):\n")
cat("  effective_drug_weeks:", round(100 * mean(estimates[, 3] < 0.05), 1), "%\n")
cat("  Interaction:", round(100 * mean(estimates[, 4] < 0.05), 1), "%\n")
