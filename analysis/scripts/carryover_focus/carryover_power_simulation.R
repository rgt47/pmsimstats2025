# Carryover-Focused Power Simulation -- Hybrid Design Only
#
# Dual-scenario comparison across two tsd_min regimes:
#   tsd_min = 0   : original (first off-drug measurement at tsd = 0)
#   tsd_min = 0.5 : corrected (minimum 0.5-week washout gap)
#
# Within each DGP, two analysis models are compared:
#   Modeled: Dbc (Hendrickson binary predictor with exponential decay)
#   Naive:   treatment (binary on/off, ignores carryover)
#
# Based on: Hendrickson et al. (2020)
# Source:   simulation_clustered.R (lines referenced in PLAN.md)

# ============================================================================
# SECTION 1: SETUP
# ============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(MASS)
  library(zoo)
  library(future)
  library(furrr)
  library(patchwork)
  library(conflicted)
})

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================================
# SECTION 2: FIXED PARAMETERS
# ============================================================================

n_participants <- 70
n_iterations <- 200

env_iterations <- Sys.getenv("N_ITERATIONS", "")
if (nchar(env_iterations) > 0) {
  n_iterations <- as.integer(env_iterations)
  cat("*** CUSTOM ITERATIONS:", n_iterations, "***\n\n")
}

biomarker_moderation <- 0.4
biomarker_correlation <- 0.3

BR_rate <- 0.5
ER_rate <- 0.2
TR_rate <- 0.1

baseline_mean <- 10.0
between_subject_sd <- 2.0
within_subject_sd <- 1.8
biomarker_mean <- 5.0
biomarker_sd <- 2.0

c.br <- 0.75
c.er <- 0.75
c.tr <- 0.75
c.cf1t <- 0.12
c.cfct <- 0.05
c.bm_baseline <- 0.25
c.baseline_resp <- 0.3

measurement_weeks_hybrid <- c(4, 8, 9, 10, 11, 12, 16, 20)

# tsd_min: minimum time-since-discontinuation (weeks)
# 0   : original Hendrickson logic. The first off-drug measurement
#       has tsd = 0, producing the (1/2)^(0/t_half) = 1 discontinuity.
# 0.5 : corrected. Imposes a 0.5-week (~3.5 day) minimum washout
#       gap, eliminating the singularity and restoring a smooth
#       power curve as a function of the carryover half-life.
tsd_min_values <- c(0, 0.5)

n_cores <- max(1, parallel::detectCores() - 1)

# ============================================================================
# SECTION 3: HELPER FUNCTIONS
# ============================================================================

build_sigma_guaranteed_pd <- function(weeks, c.bm, params) {
  n_tp <- length(weeks)

  br_idx <- 1:n_tp
  er_idx <- (n_tp + 1):(2 * n_tp)
  tr_idx <- (2 * n_tp + 1):(3 * n_tp)
  bm_idx <- 3 * n_tp + 1
  bl_idx <- 3 * n_tp + 2

  Sigma_22 <- matrix(c(
    biomarker_sd^2,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    between_subject_sd^2
  ), 2, 2)

  build_ar1_time <- function(weeks, rho, sigma) {
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

  Sigma_BR <- build_ar1_time(weeks, c.br, within_subject_sd)
  Sigma_ER <- build_ar1_time(weeks, c.er, within_subject_sd)
  Sigma_TR <- build_ar1_time(weeks, c.tr, within_subject_sd)

  Sigma_11 <- matrix(0, 3 * n_tp, 3 * n_tp)
  Sigma_11[br_idx, br_idx] <- Sigma_BR
  Sigma_11[er_idx, er_idx] <- Sigma_ER
  Sigma_11[tr_idx, tr_idx] <- Sigma_TR

  for (i in 1:n_tp) {
    for (j in 1:n_tp) {
      time_lag <- abs(weeks[i] - weeks[j])
      cross_cov <- if (i == j) {
        c.cf1t * within_subject_sd^2
      } else {
        c.cfct * within_subject_sd^2 * (0.9^time_lag)
      }
      Sigma_11[br_idx[i], er_idx[j]] <- cross_cov
      Sigma_11[er_idx[j], br_idx[i]] <- cross_cov
      Sigma_11[br_idx[i], tr_idx[j]] <- cross_cov
      Sigma_11[tr_idx[j], br_idx[i]] <- cross_cov
      Sigma_11[er_idx[i], tr_idx[j]] <- cross_cov
      Sigma_11[tr_idx[j], er_idx[i]] <- cross_cov
    }
  }

  effective_c.bm <- c.bm
  effective_c.baseline <- c.baseline_resp

  Sigma_12 <- matrix(0, 3 * n_tp, 2)
  Sigma_12[br_idx, 1] <- effective_c.bm * within_subject_sd * biomarker_sd
  Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd *
    biomarker_sd
  Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd *
    biomarker_sd
  Sigma_12[br_idx, 2] <- effective_c.baseline * within_subject_sd *
    between_subject_sd
  Sigma_12[er_idx, 2] <- effective_c.baseline * within_subject_sd *
    between_subject_sd
  Sigma_12[tr_idx, 2] <- effective_c.baseline * within_subject_sd *
    between_subject_sd

  Sigma_22_inv <- solve(Sigma_22)
  cross_term <- Sigma_12 %*% Sigma_22_inv %*% t(Sigma_12)
  Sigma_cond <- Sigma_11 - cross_term
  min_eig <- min(eigen(Sigma_cond, only.values = TRUE)$values)

  if (min_eig <= 1e-6) {
    stop(sprintf(
      "Sigma not positive definite (min eigenvalue = %.2e)", min_eig
    ))
  }

  total_dim <- 3 * n_tp + 2
  Sigma <- matrix(0, total_dim, total_dim)
  Sigma[1:(3 * n_tp), 1:(3 * n_tp)] <- Sigma_11
  Sigma[(3 * n_tp + 1):total_dim, (3 * n_tp + 1):total_dim] <- Sigma_22
  Sigma[1:(3 * n_tp), (3 * n_tp + 1):total_dim] <- Sigma_12
  Sigma[(3 * n_tp + 1):total_dim, 1:(3 * n_tp)] <- t(Sigma_12)

  indices <- list(
    br = br_idx, er = er_idx, tr = tr_idx,
    bm = bm_idx, bl = bl_idx, n_tp = n_tp
  )

  cond_mean_transform <- Sigma_12 %*% Sigma_22_inv

  list(
    Sigma = Sigma,
    Sigma_22 = Sigma_22,
    Sigma_22_inv = Sigma_22_inv,
    Sigma_cond = Sigma_cond,
    Sigma_12 = Sigma_12,
    cond_mean_transform = cond_mean_transform,
    indices = indices,
    effective_c.bm = effective_c.bm
  )
}

generate_participant_twostage <- function(sigma_parts, idx) {
  Sigma_22 <- sigma_parts$Sigma_22
  Sigma_cond <- sigma_parts$Sigma_cond
  n_tp <- idx$n_tp

  stage1 <- mvrnorm(
    1, mu = c(biomarker_mean, baseline_mean), Sigma = Sigma_22
  )
  biomarker <- stage1[1]
  baseline <- stage1[2]

  z <- c(biomarker - biomarker_mean, baseline - baseline_mean)
  cond_mean <- as.vector(sigma_parts$cond_mean_transform %*% z)
  responses <- mvrnorm(1, mu = cond_mean, Sigma = Sigma_cond)

  list(
    biomarker = biomarker,
    baseline = baseline,
    br_random = responses[idx$br],
    er_random = responses[idx$er],
    tr_random = responses[idx$tr]
  )
}

create_hybrid_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))

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
      ),
      expectancy = if_else(week %in% c(4, 8), 1.0, 0.5)
    )
}

# ============================================================================
# SECTION 4: PARAMETER SWEEP
# ============================================================================
# Three dimensions:
#   - tsd_min: 0 (uncorrected) vs 0.5 (corrected washout gap)
#   - carryover_halflife: 5 values (finer grid to reveal the shape)
#   - model_carryover: TRUE (Dbc) vs FALSE (treatment)
# At carryover=0 both predictors are identical, so we keep only one
# row there per tsd_min value.

param_grid <- expand_grid(
  tsd_min = tsd_min_values,
  carryover_halflife = c(0, 0.05, 0.1, 0.15, 0.2),
  model_carryover = c(TRUE, FALSE)
) |>
  filter(!(carryover_halflife == 0 & model_carryover == FALSE)) |>
  mutate(condition_id = seq_len(n()))

cat(sprintf(
  "Carryover-focused simulation: %d conditions x %d iterations\n",
  nrow(param_grid), n_iterations
))
cat(sprintf(
  "Design: Hybrid | n=%d | bm_mod=%.1f | bm_corr=%.1f\n",
  n_participants, biomarker_moderation, biomarker_correlation
))
cat(sprintf(
  "tsd_min values: %s\n",
  paste(tsd_min_values, collapse = " vs ")
))
cat("Modeled:  Dbc * bm_centered + week\n")
cat("Naive:    treatment * bm_centered + week\n\n")

# ============================================================================
# SECTION 5: SINGLE ITERATION FUNCTION
# ============================================================================

sigma_parts <- build_sigma_guaranteed_pd(
  weeks = measurement_weeks_hybrid,
  c.bm = biomarker_correlation,
  params = list()
)
idx <- sigma_parts$indices

cat("Sigma matrix validated (PD check passed)\n\n")

run_single_iteration <- function(iter, params) {
  set.seed(iter * 1000 + params$condition_id)

  local_n <- n_participants
  bm_mod <- biomarker_moderation
  carryover_t_half <- params$carryover_halflife
  local_tsd_min <- params$tsd_min

  trial_design <- create_hybrid_design(local_n, measurement_weeks_hybrid)

  all_participant_data <- lapply(1:local_n, function(pid) {
    result <- generate_participant_twostage(sigma_parts, idx)
    tibble(
      participant_id = pid,
      biomarker = result$biomarker,
      baseline = result$baseline,
      br_random = result$br_random,
      er_random = result$er_random,
      tr_random = result$tr_random
    )
  })

  participant_data <- bind_rows(all_participant_data) |>
    group_by(participant_id) |>
    mutate(timepoint_idx = row_number()) |>
    ungroup()

  trial_data <- trial_design |>
    group_by(participant_id) |>
    mutate(timepoint_idx = row_number()) |>
    ungroup() |>
    left_join(participant_data,
              by = c("participant_id", "timepoint_idx")) |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      weeks_with_expectancy = cumsum(expectancy),
      weeks_in_trial = week - min(week),
      just_discontinued = treatment == 0 &
        lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(
        just_discontinued, week, NA_real_
      ),
      discontinuation_week = zoo::na.locf(
        discontinuation_week, na.rm = FALSE
      ),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        pmax(week - discontinuation_week, local_tsd_min),
        0
      ),
      bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
      effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered),
      br_on_drug = weeks_on_drug * effective_BR_rate,
      br_at_discont = if_else(
        just_discontinued, lag(br_on_drug, default = 0), NA_real_
      ),
      br_at_discont = zoo::na.locf(br_at_discont, na.rm = FALSE),
      Dbc = case_when(
        treatment == 1 ~ 1,
        treatment == 0 & carryover_t_half == 0 ~ 0,
        treatment == 0 & carryover_t_half > 0 ~
          (1 / 2)^(tsd / carryover_t_half),
        TRUE ~ 0
      ),
      br_off_drug = if_else(
        treatment == 0 & !is.na(br_at_discont),
        br_at_discont * Dbc,
        0
      ),
      BR_mean = if_else(treatment == 1, br_on_drug, br_off_drug),
      BR = BR_mean + br_random,
      ER_mean = weeks_with_expectancy * ER_rate,
      ER = ER_mean + er_random,
      TR_mean = weeks_in_trial * TR_rate,
      TR = TR_mean + tr_random,
      response = baseline + BR + ER + TR
    ) |>
    ungroup()

  trial_data <- trial_data |>
    mutate(
      bm_centered = (biomarker - mean(biomarker)) / sd(biomarker)
    )

  model_result <- tryCatch({
    if (params$model_carryover) {
      model <- lme(
        response ~ Dbc * bm_centered + week,
        random = ~ 1 | participant_id,
        correlation = corCAR1(form = ~ week | participant_id),
        data = trial_data,
        control = lmeControl(opt = "optim")
      )
      interaction_term <- "Dbc:bm_centered"
    } else {
      model <- lme(
        response ~ treatment * bm_centered + week,
        random = ~ 1 | participant_id,
        correlation = corCAR1(form = ~ week | participant_id),
        data = trial_data,
        control = lmeControl(opt = "optim")
      )
      interaction_term <- "treatment:bm_centered"
    }

    coefs <- summary(model)$tTable
    idx_coef <- which(rownames(coefs) == interaction_term)

    tibble(
      iteration = iter,
      tsd_min = local_tsd_min,
      carryover_halflife = params$carryover_halflife,
      model_carryover = params$model_carryover,
      effect_size = coefs[idx_coef, "Value"],
      se = coefs[idx_coef, "Std.Error"],
      t_value = coefs[idx_coef, "t-value"],
      p_value = coefs[idx_coef, "p-value"],
      significant = coefs[idx_coef, "p-value"] < 0.05
    )
  }, error = function(e) {
    tibble(
      iteration = iter,
      tsd_min = local_tsd_min,
      carryover_halflife = params$carryover_halflife,
      model_carryover = params$model_carryover,
      effect_size = NA_real_,
      se = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      significant = NA
    )
  })

  model_result
}

# ============================================================================
# SECTION 6: SIMULATION LOOP
# ============================================================================

plan(multisession, workers = n_cores)
cat(sprintf("Parallel backend: %d workers\n\n", n_cores))

start_time <- Sys.time()
results_list <- list()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf(
    "Condition %d/%d: tsd_min=%.1f, co_t1/2=%.2f wks, model=%s\n",
    i, nrow(param_grid),
    params$tsd_min,
    params$carryover_halflife,
    if (params$model_carryover) "modeled" else "naive"
  ))

  local_params <- params
  local_sigma_parts <- sigma_parts
  local_idx <- idx

  condition_results <- future_map_dfr(
    1:n_iterations,
    function(iter) {
      run_single_iteration(iter, local_params)
    },
    .options = furrr_options(
      seed = TRUE,
      globals = list(
        run_single_iteration = run_single_iteration,
        generate_participant_twostage = generate_participant_twostage,
        create_hybrid_design = create_hybrid_design,
        build_sigma_guaranteed_pd = build_sigma_guaranteed_pd,
        sigma_parts = local_sigma_parts,
        idx = local_idx,
        local_params = local_params,
        n_participants = n_participants,
        biomarker_moderation = biomarker_moderation,
        biomarker_correlation = biomarker_correlation,
        biomarker_mean = biomarker_mean,
        biomarker_sd = biomarker_sd,
        baseline_mean = baseline_mean,
        between_subject_sd = between_subject_sd,
        within_subject_sd = within_subject_sd,
        measurement_weeks_hybrid = measurement_weeks_hybrid,
        BR_rate = BR_rate,
        ER_rate = ER_rate,
        TR_rate = TR_rate,
        c.br = c.br,
        c.er = c.er,
        c.tr = c.tr,
        c.cf1t = c.cf1t,
        c.cfct = c.cfct,
        c.bm_baseline = c.bm_baseline,
        c.baseline_resp = c.baseline_resp
      ),
      packages = c("tidyverse", "nlme", "MASS", "zoo")
    )
  )

  results_list[[i]] <- condition_results

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  completed <- i * n_iterations
  total <- nrow(param_grid) * n_iterations
  rate <- completed / elapsed
  remaining <- (total - completed) / rate
  cat(sprintf(
    "  Done (%.0fs elapsed, ~%.0fs remaining)\n", elapsed, remaining
  ))
}

results <- bind_rows(results_list)
plan(sequential)

# ============================================================================
# SECTION 7: RESULTS SUMMARY
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("RESULTS: Carryover Impact on Power (Hybrid Design)\n")
cat(strrep("=", 70), "\n\n")

summary_results <- results |>
  group_by(tsd_min, carryover_halflife, model_carryover) |>
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    n = n(),
    n_errors = sum(is.na(significant)),
    .groups = "drop"
  ) |>
  mutate(
    ci_lower = pmax(0,
      (power + 1.96^2 / (2 * n) -
        1.96 * sqrt(
          (power * (1 - power) + 1.96^2 / (4 * n)) / n
        )) / (1 + 1.96^2 / n)
    ),
    ci_upper = pmin(1,
      (power + 1.96^2 / (2 * n) +
        1.96 * sqrt(
          (power * (1 - power) + 1.96^2 / (4 * n)) / n
        )) / (1 + 1.96^2 / n)
    ),
    scenario = if_else(
      model_carryover, "Modeled (Dbc)", "Naive (treatment)"
    )
  )

for (tv in tsd_min_values) {
  label <- if (tv == 0) "UNCORRECTED (tsd_min = 0)"
    else sprintf("CORRECTED (tsd_min = %.1f)", tv)
  cat(sprintf("\n--- %s ---\n", label))
  print(
    summary_results |>
      filter(tsd_min == tv) |>
      select(carryover_halflife, scenario, power, ci_lower,
             ci_upper, mean_effect, sd_effect, n_errors),
    n = 20
  )
}

scenario_comparison <- summary_results |>
  filter(carryover_halflife > 0) |>
  select(tsd_min, carryover_halflife,
         model_carryover, power) |>
  pivot_wider(
    names_from = model_carryover,
    values_from = power,
    names_prefix = "model_"
  ) |>
  mutate(
    power_difference = model_TRUE - model_FALSE,
    pct_change = sprintf(
      "%+.0f%%", (power_difference / model_FALSE) * 100
    )
  )

for (tv in tsd_min_values) {
  label <- if (tv == 0) "UNCORRECTED (tsd_min = 0)"
    else sprintf("CORRECTED (tsd_min = %.1f)", tv)
  cat(sprintf("\nModeled vs naive comparison -- %s:\n", label))
  print(
    scenario_comparison |> filter(tsd_min == tv),
    n = 10
  )
}

# ============================================================================
# SECTION 8: VISUALIZATION
# ============================================================================

tsd_labels <- c(
  "0" = "Uncorrected (tsd_min = 0)",
  "0.5" = "Corrected (tsd_min = 0.5 wks)"
)

plot_data <- summary_results |>
  mutate(
    scenario = factor(scenario,
      levels = c("Modeled (Dbc)", "Naive (treatment)")
    ),
    tsd_label = factor(
      tsd_labels[as.character(tsd_min)],
      levels = tsd_labels
    )
  )

provenance <- sprintf(
  "Source: carryover_power_simulation.R | %s | n=%d, iter=%d, bm_mod=%.1f",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  n_participants, n_iterations, biomarker_moderation
)

p_power <- ggplot(
  plot_data,
  aes(x = carryover_halflife, y = power,
      color = scenario, shape = scenario)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.008, linewidth = 0.6
  ) +
  geom_hline(yintercept = 0.80, linetype = "dashed", alpha = 0.4) +
  facet_wrap(~ tsd_label) +
  scale_color_manual(values = c(
    "Modeled (Dbc)" = "#2ca02c",
    "Naive (treatment)" = "#d62728"
  )) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Statistical Power",
    x = "Carryover half-life (weeks)",
    y = "Power (proportion significant)",
    color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_effect <- ggplot(
  plot_data,
  aes(x = carryover_halflife, y = mean_effect,
      color = scenario, shape = scenario)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_effect - sd_effect,
        ymax = mean_effect + sd_effect),
    width = 0.008, linewidth = 0.6, alpha = 0.5
  ) +
  facet_wrap(~ tsd_label) +
  scale_color_manual(values = c(
    "Modeled (Dbc)" = "#2ca02c",
    "Naive (treatment)" = "#d62728"
  )) +
  labs(
    title = "Mean Effect Size",
    x = "Carryover half-life (weeks)",
    y = "Mean interaction estimate",
    color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_combined <- p_power / p_effect +
  plot_annotation(
    title = "Carryover Effect on Biomarker Moderation Detection",
    subtitle = "Hybrid Design, N=70, biomarker moderation=0.4",
    caption = provenance,
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 8, color = "grey50")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "carryover_power_hybrid.pdf"),
  p_combined, width = 11, height = 10
)
ggsave(
  file.path(output_dir, "carryover_power_hybrid.png"),
  p_combined, width = 11, height = 10, dpi = 300
)

cat(sprintf("\nFigure saved to %s\n", output_dir))

# ============================================================================
# SECTION 9: SAVE RESULTS
# ============================================================================

sim_parameters <- list(
  n_participants = n_participants,
  n_iterations = n_iterations,
  biomarker_moderation = biomarker_moderation,
  biomarker_correlation = biomarker_correlation,
  design = "hybrid",
  carryover_halflife_values = c(0, 0.05, 0.1, 0.15, 0.2),
  tsd_min_values = tsd_min_values,
  model_modeled = "Dbc * bm_centered + week",
  model_naive = "treatment * bm_centered + week"
)

save(
  results, summary_results, scenario_comparison, sim_parameters,
  file = file.path(output_dir, "carryover_focus_results.RData")
)

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("Output files:\n")
cat("  - carryover_power_hybrid.pdf/png (4-panel: uncorrected vs corrected)\n")
cat("  - carryover_focus_results.RData\n")
cat(sprintf("Total time: %.1f minutes\n",
  as.numeric(difftime(Sys.time(), start_time, units = "mins"))
))
