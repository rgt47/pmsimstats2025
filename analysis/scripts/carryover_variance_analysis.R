# Analysis: How Carryover Reduces Variance in effective_drug_weeks
# This explains why power drops with longer carryover half-lives

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

set.seed(42)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("HOW CARRYOVER REDUCES VARIANCE IN effective_drug_weeks\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Hybrid design measurement weeks
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)
n_participants <- 70

# Create one example of hybrid design
path_assignment <- sample(rep(1:4, length.out = n_participants))

create_design_data <- function(carryover_t_half) {
  expand_grid(
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
      )
    ) |>
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
      weeks_at_discont = if_else(just_discontinued, lag(weeks_on_drug, default = 0), NA_real_),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE),
      carryover_effect = if_else(
        treatment == 0 & carryover_t_half > 0,
        (1/2)^(tsd / carryover_t_half),
        0
      ),
      effective_drug_weeks = if_else(
        treatment == 1,
        weeks_on_drug,
        if_else(!is.na(weeks_at_discont) & carryover_t_half > 0,
                weeks_at_discont * carryover_effect,
                0)
      )
    ) |>
    ungroup() |>
    mutate(carryover_t_half = carryover_t_half)
}

# Compare across carryover levels
carryover_levels <- c(0, 0.1, 0.2)  # Hendrickson values (weeks)
all_data <- map_dfr(carryover_levels, create_design_data)

cat("VARIANCE OF effective_drug_weeks BY CARRYOVER HALF-LIFE\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

variance_summary <- all_data |>
  group_by(carryover_t_half) |>
  summarise(
    mean_edw = mean(effective_drug_weeks),
    var_edw = var(effective_drug_weeks),
    sd_edw = sd(effective_drug_weeks),
    min_edw = min(effective_drug_weeks),
    max_edw = max(effective_drug_weeks),
    range_edw = max_edw - min_edw,
    .groups = "drop"
  )

print(variance_summary)

cat("\n")
cat("EXPLANATION:\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("With NO carryover (t_half = 0):\n")
cat("  - On drug: effective_drug_weeks = weeks_on_drug (1, 2, 3, 4, 5...)\n")
cat("  - Off drug: effective_drug_weeks = 0\n")
cat("  - This creates MAXIMUM contrast between on/off drug states\n\n")

cat("With carryover (t_half > 0):\n")
cat("  - On drug: effective_drug_weeks = weeks_on_drug (same as before)\n")
cat("  - Off drug: effective_drug_weeks = weeks_at_discont * decay > 0\n")
cat("  - The off-drug values are NO LONGER ZERO!\n")
cat("  - This REDUCES the contrast between on/off drug states\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DETAILED EXAMPLE: ONE PARTICIPANT (Path 3)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

example_participant <- all_data |>
  filter(participant_id == 1, path == 3) |>
  dplyr::select(week, treatment, weeks_on_drug, tsd, carryover_effect,
         effective_drug_weeks, carryover_t_half) |>
  pivot_wider(
    id_cols = c(week, treatment, weeks_on_drug),
    names_from = carryover_t_half,
    values_from = effective_drug_weeks,
    names_prefix = "edw_t"
  )

# Find a path 3 participant
path3_id <- which(path_assignment == 3)[1]

example_data <- all_data |>
  filter(participant_id == path3_id) |>
  dplyr::select(week, treatment, weeks_on_drug, tsd, carryover_t_half,
         carryover_effect, effective_drug_weeks)

cat("Path 3 treatment schedule:\n")
cat("  Week 4-9: On drug\n")
cat("  Week 10: Off drug (discontinue)\n")
cat("  Week 11-12: Off drug\n")
cat("  Week 16: On drug (re-start)\n")
cat("  Week 20: Off drug\n\n")

for (t_half in carryover_levels) {
  cat(sprintf("Carryover t_half = %.1f weeks:\n", t_half))
  example_data |>
    filter(carryover_t_half == t_half) |>
    dplyr::select(week, treatment, weeks_on_drug, tsd, carryover_effect,
           effective_drug_weeks) |>
    print(n = 8)
  cat("\n")
}

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("WHY REDUCED VARIANCE HURTS POWER\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("In regression, the standard error of a coefficient is:\n\n")
cat("  SE(β) = σ / sqrt(Σ(x_i - x̄)²)\n\n")
cat("Where the denominator is the SUM OF SQUARED DEVIATIONS of the predictor.\n\n")

cat("When variance of effective_drug_weeks decreases:\n")
cat("  - Σ(x_i - x̄)² decreases\n")
cat("  - SE(β) INCREASES\n")
cat("  - t-statistic = β/SE(β) DECREASES\n")
cat("  - Power to detect the effect DECREASES\n\n")

cat("This is a fundamental statistical principle:\n")
cat("  Less variation in the predictor → Less precision in estimating its effect\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("THE KEY INSIGHT\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("With no carryover:\n")
cat("  - Off-drug observations have effective_drug_weeks = 0\n")
cat("  - On-drug observations have effective_drug_weeks = 1, 2, 3, 4, 5...\n")
cat("  - BIG difference between on and off states\n\n")

cat("With long carryover (e.g., t_half = 2 weeks):\n")
cat("  - Off-drug observations have effective_drug_weeks = 3*0.84, 3*0.71...\n")
cat("  - These values are CLOSE to the on-drug values\n")
cat("  - SMALL difference between on and off states\n\n")

cat("The carryover effect 'fills in' the zero values during off-drug periods,\n")
cat("making the data look more like a monotonic trend rather than a contrast\n")
cat("between distinct treatment states.\n\n")

# Visualize
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CORRELATION WITH WEEK (COLLINEARITY)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cor_summary <- all_data |>
  group_by(carryover_t_half) |>
  summarise(
    cor_with_week = cor(effective_drug_weeks, week),
    .groups = "drop"
  )

print(cor_summary)

cat("\nAs carryover increases, effective_drug_weeks becomes MORE correlated with\n")
cat("week (time). This collinearity further reduces the unique information\n")
cat("available to estimate the treatment effect.\n")
