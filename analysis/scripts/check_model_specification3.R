# Model Specification Diagnostic - Part 3
# Compare lmer (no AR1) vs lme (with AR1) and check bm_centered scaling

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(nlme)
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
between_subject_sd <- 2.0
within_subject_sd <- 1.8
bm_mod <- 0.65

c.br <- 0.75
c.er <- 0.75
c.tr <- 0.75

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL SPECIFICATION COMPARISON: lmer vs lme\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)
carryover_t_half <- 1.0

# Build AR(1) covariance for random effects (simplified version)
build_ar1_cov <- function(weeks, rho, sigma) {
  n <- length(weeks)
  Cov <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      time_lag <- abs(weeks[i] - weeks[j])
      Cov[i, j] <- sigma^2 * rho^time_lag
    }
  }
  Cov
}

generate_participant <- function(pid, weeks, bm_mod, carryover_t_half) {
  n_tp <- length(weeks)

  biomarker <- rnorm(1, biomarker_mean, biomarker_sd)
  baseline <- rnorm(1, baseline_mean, between_subject_sd)
  random_intercept <- rnorm(1, 0, between_subject_sd)

  Sigma_BR <- build_ar1_cov(weeks, c.br, within_subject_sd)
  Sigma_ER <- build_ar1_cov(weeks, c.er, within_subject_sd)
  Sigma_TR <- build_ar1_cov(weeks, c.tr, within_subject_sd)

  br_random <- mvrnorm(1, rep(0, n_tp), Sigma_BR)
  er_random <- mvrnorm(1, rep(0, n_tp), Sigma_ER)
  tr_random <- mvrnorm(1, rep(0, n_tp), Sigma_TR)

  tibble(
    participant_id = pid,
    biomarker = biomarker,
    baseline = baseline,
    random_intercept = random_intercept,
    timepoint_idx = 1:n_tp,
    br_random = br_random,
    er_random = er_random,
    tr_random = tr_random
  )
}

path_assignment <- sample(rep(1:4, length.out = n_participants))

trial_design <- expand_grid(
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
  ) |>
  group_by(participant_id) |>
  mutate(timepoint_idx = row_number()) |>
  ungroup()

participant_data <- map_dfr(1:n_participants, ~generate_participant(
  .x, measurement_weeks, bm_mod, carryover_t_half
))

trial_data <- trial_design |>
  left_join(participant_data, by = c("participant_id", "timepoint_idx")) |>
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
    BR = BR_mean + br_random,
    ER_mean = weeks_with_expectancy * ER_rate,
    ER = ER_mean + er_random,
    TR_mean = weeks_in_trial * TR_rate,
    TR = TR_mean + tr_random,
    response = baseline + random_intercept + BR + ER + TR
  ) |>
  ungroup()

# Create different versions of bm_centered
trial_data <- trial_data |>
  mutate(
    bm_centered_sample = biomarker - mean(biomarker),
    bm_centered_std = (biomarker - mean(biomarker)) / sd(biomarker)
  )

cat("TRUE PARAMETERS:\n")
cat("  BR_rate:", BR_rate, "\n")
cat("  bm_mod:", bm_mod, "\n")
cat("  BR_rate * bm_mod:", BR_rate * bm_mod, "(standardized scale)\n")
cat("  BR_rate * bm_mod / biomarker_sd:", BR_rate * bm_mod / biomarker_sd,
    "(unstandardized scale)\n\n")

cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON 1: lmer vs lme (with standardized bm_centered)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# lmer (no AR1)
model_lmer <- lmer(
  response ~ effective_drug_weeks * bm_centered_std + week + (1 | participant_id),
  data = trial_data,
  REML = TRUE
)

# lme with AR1
model_lme <- lme(
  response ~ effective_drug_weeks * bm_centered_std + week,
  random = ~ 1 | participant_id,
  correlation = corCAR1(form = ~ week | participant_id),
  data = trial_data,
  control = lmeControl(opt = "optim")
)

cat("lmer (no AR1 correlation):\n")
print(summary(model_lmer)$coefficients)
cat("\n")

cat("lme (with corCAR1):\n")
print(summary(model_lme)$tTable)
cat("\n")

# Extract AR1 correlation estimate
phi <- coef(model_lme$modelStruct$corStruct, unconstrained = FALSE)
cat("Estimated AR1 correlation (phi):", round(phi, 4), "\n")
cat("True AR1 correlation (c.br, c.er, c.tr):", c.br, "\n\n")

cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON 2: Effect of bm_centered scaling\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# Using unstandardized (current simulation approach)
model_lme_unstd <- lme(
  response ~ effective_drug_weeks * bm_centered_sample + week,
  random = ~ 1 | participant_id,
  correlation = corCAR1(form = ~ week | participant_id),
  data = trial_data,
  control = lmeControl(opt = "optim")
)

cat("lme with UNSTANDARDIZED bm_centered (current simulation):\n")
cat("  Expected interaction coef: BR_rate * bm_mod / biomarker_sd =",
    BR_rate * bm_mod / biomarker_sd, "\n")
cat("  Estimated:", round(fixef(model_lme_unstd)["effective_drug_weeks:bm_centered_sample"], 4), "\n\n")

cat("lme with STANDARDIZED bm_centered:\n")
cat("  Expected interaction coef: BR_rate * bm_mod =", BR_rate * bm_mod, "\n")
cat("  Estimated:", round(fixef(model_lme)["effective_drug_weeks:bm_centered_std"], 4), "\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SIMULATION: COMPARE POWER WITH CORRECT SPECIFICATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

n_sim <- 50
results <- tibble(
  iteration = integer(),
  model_type = character(),
  estimate = double(),
  se = double(),
  p_value = double()
)

cat("Running", n_sim, "iterations...\n")

for (i in 1:n_sim) {
  set.seed(i * 1000)

  path_assignment <- sample(rep(1:4, length.out = n_participants))

  trial_design_i <- expand_grid(
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
    ) |>
    group_by(participant_id) |>
    mutate(timepoint_idx = row_number()) |>
    ungroup()

  participant_data_i <- map_dfr(1:n_participants, ~generate_participant(
    .x, measurement_weeks, bm_mod, carryover_t_half
  ))

  trial_data_i <- trial_design_i |>
    left_join(participant_data_i, by = c("participant_id", "timepoint_idx")) |>
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
      BR = BR_mean + br_random,
      ER_mean = weeks_with_expectancy * ER_rate,
      ER = ER_mean + er_random,
      TR_mean = weeks_in_trial * TR_rate,
      TR = TR_mean + tr_random,
      response = baseline + random_intercept + BR + ER + TR
    ) |>
    ungroup() |>
    mutate(bm_centered_std = (biomarker - mean(biomarker)) / sd(biomarker))

  # lmer
  m_lmer <- tryCatch({
    lmer(
      response ~ effective_drug_weeks * bm_centered_std + week + (1 | participant_id),
      data = trial_data_i,
      REML = TRUE
    )
  }, error = function(e) NULL)

  if (!is.null(m_lmer)) {
    coefs <- summary(m_lmer)$coefficients
    t_val <- coefs["effective_drug_weeks:bm_centered_std", "t value"]
    results <- bind_rows(results, tibble(
      iteration = i,
      model_type = "lmer",
      estimate = coefs["effective_drug_weeks:bm_centered_std", "Estimate"],
      se = coefs["effective_drug_weeks:bm_centered_std", "Std. Error"],
      p_value = 2 * pnorm(-abs(t_val))
    ))
  }

  # lme with AR1
  m_lme <- tryCatch({
    lme(
      response ~ effective_drug_weeks * bm_centered_std + week,
      random = ~ 1 | participant_id,
      correlation = corCAR1(form = ~ week | participant_id),
      data = trial_data_i,
      control = lmeControl(opt = "optim")
    )
  }, error = function(e) NULL)

  if (!is.null(m_lme)) {
    coefs <- summary(m_lme)$tTable
    results <- bind_rows(results, tibble(
      iteration = i,
      model_type = "lme_AR1",
      estimate = coefs["effective_drug_weeks:bm_centered_std", "Value"],
      se = coefs["effective_drug_weeks:bm_centered_std", "Std.Error"],
      p_value = coefs["effective_drug_weeks:bm_centered_std", "p-value"]
    ))
  }
}

cat("\nResults summary:\n")
results |>
  group_by(model_type) |>
  summarise(
    n = n(),
    mean_estimate = mean(estimate),
    sd_estimate = sd(estimate),
    mean_se = mean(se),
    power = mean(p_value < 0.05),
    .groups = "drop"
  ) |>
  print()

cat("\nTrue interaction coefficient: BR_rate * bm_mod =", BR_rate * bm_mod, "\n")
cat("\nConclusion:\n")
cat("- If estimates are unbiased, mean should be near", BR_rate * bm_mod, "\n")
cat("- lmer SE may be incorrect due to ignoring AR1 correlation\n")
cat("- lme with AR1 should give more accurate SE estimates\n")
