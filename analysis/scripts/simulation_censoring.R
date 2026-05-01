# ============================================================================
# CENSORING SENSITIVITY ANALYSIS
# ============================================================================
# Replicates Hendrickson Figure 4: Power as a function of trial design,
# biomarker effect, carryover, and censoring/dropout patterns
#
# Censoring model (Hendrickson p.6):
#   P(dropout per timepoint) = β₀ + β₁ * (change_from_baseline + 100)²
#   - β₀: flat hazard (baseline dropout risk)
#   - β₁: symptom-dependent dropout (higher for worse symptoms)
#
# Usage: Rscript simulation_censoring.R
#        EXPLORATORY_MODE=TRUE Rscript simulation_censoring.R
#        N_ITERATIONS=100 Rscript simulation_censoring.R
# ============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(lme4)
  library(MASS)
  library(zoo)
  library(patchwork)
  library(conflicted)
  library(future)
  library(furrr)
})

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(nlme::lme)

# ============================================================================
# PARAMETERS
# ============================================================================

n_iterations <- 500

# Environment variable overrides
exploratory_mode <- as.logical(Sys.getenv("EXPLORATORY_MODE", "FALSE"))
exploratory_iterations <- 20

env_iterations <- Sys.getenv("N_ITERATIONS", "")
if (nchar(env_iterations) > 0) {
  n_iterations <- as.integer(env_iterations)
  cat("*** CUSTOM ITERATIONS:", n_iterations, "***\n\n")
} else if (exploratory_mode) {
  n_iterations <- exploratory_iterations
  cat("*** EXPLORATORY MODE: Running with", n_iterations, "iterations ***\n\n")
}

# Parallel processing
n_cores <- max(1, parallel::detectCores() - 1)
use_parallel <- TRUE
model_engine <- "lme"

# Three-factor response model
BR_rate <- 0.5
ER_rate <- 0.2
TR_rate <- 0.1
treatment_effect <- BR_rate

baseline_mean <- 10.0
between_subject_sd <- 2.0
within_subject_sd <- 1.8
biomarker_mean <- 5.0
biomarker_sd <- 2.0

# Correlation parameters
c.br <- 0.75
c.er <- 0.75
c.tr <- 0.75
c.cf1t <- 0.12
c.cfct <- 0.05
c.bm_baseline <- 0.25
c.baseline_resp <- 0.3

carryover_mode <- "halflife"

# ============================================================================
# CENSORING PARAMETERS (Hendrickson-style, rescaled for our response scale)
# ============================================================================
# P(dropout) = β₀ + β₁ * (symptom_change + offset)²
# where symptom_change = response - baseline
#
# Hendrickson used offset=100 for PTSD scale (0-136). Our response scale is
# much smaller (~10-30), so we use offset=20 and much smaller β₁ values.
#
# Target dropout rates (cumulative over 8 timepoints):
#   - none: 0%
#   - balanced: ~15-25%
#   - more_flat: ~30-40%
#   - more_biased: ~20-30% (biased toward non-responders)
#   - high_dropout: ~40-50%

censoring_params <- tibble::tribble(
  ~censoring_type, ~beta_0,  ~beta_1,    ~description,
  "none",          0.00,     0.00,       "No censoring",
  "balanced",      0.02,     0.00002,    "Balanced dropout (~20%)",
  "more_flat",     0.04,     0.00001,    "More flat (higher baseline, ~35%)",
  "more_biased",   0.01,     0.00005,    "More biased (symptom-dependent, ~25%)",
  "high_dropout",  0.05,     0.00005,    "High dropout (~45%)"
)

# Censoring function: returns TRUE if participant drops out at this timepoint
# offset=20 for our response scale (baseline ~10, changes ±10)
apply_censoring <- function(response, baseline, beta_0, beta_1, offset = 20) {
  if (beta_0 == 0 && beta_1 == 0) {
    FALSE
  } else {
    symptom_change <- response - baseline
    p_dropout <- beta_0 + beta_1 * (symptom_change + offset)^2
    p_dropout <- pmin(pmax(p_dropout, 0), 1)
    runif(1) < p_dropout
  }
}

# ============================================================================
# CARRYOVER FUNCTION
# ============================================================================

calculate_carryover <- function(tsd, carryover_param, mode = carryover_mode) {
  if (carryover_param == 0 || is.na(tsd) || tsd <= 0) {
    0
  } else if (mode == "halflife") {
    (1/2)^(tsd / carryover_param)
  } else {
    stop("Unknown carryover mode: ", mode)
  }
}

# ============================================================================
# PARAMETER GRID - Censoring Sensitivity Analysis
# ============================================================================
# Matches Hendrickson Figure 4 structure:
# - 4 trial designs
# - 3 biomarker effect levels (0, 0.3, 0.6)
# - 3 carryover levels (0, 0.1, 0.2 weeks - using Hendrickson's shorter values)
# - 5 censoring levels
# - 2 sample sizes (35, 70)

param_grid <- expand_grid(
  design = c("ol", "ol_bdc", "crossover", "hybrid"),
  n_participants = c(35, 70),
  biomarker_moderation = c(0, 0.3, 0.6),
  carryover_halflife = c(0, 0.1, 0.2),  # Hendrickson's shorter half-lives
  censoring_type = c("none", "balanced", "more_flat", "more_biased", "high_dropout")
) |>
  left_join(censoring_params, by = "censoring_type") |>
  mutate(
    biomarker_correlation = if_else(biomarker_moderation == 0, 0, 0.3),
    condition_id = row_number()
  )

cat("Total conditions:", nrow(param_grid), "\n")
cat("Conditions breakdown:\n")
cat("  Designs:", length(unique(param_grid$design)), "\n")
cat("  Sample sizes:", paste(unique(param_grid$n_participants), collapse = ", "), "\n")
cat("  Biomarker moderation:", paste(unique(param_grid$biomarker_moderation), collapse = ", "), "\n")
cat("  Carryover half-lives:", paste(unique(param_grid$carryover_halflife), collapse = ", "), "\n")
cat("  Censoring types:", paste(unique(param_grid$censoring_type), collapse = ", "), "\n\n")

# ============================================================================
# SIGMA MATRIX BUILDER
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
  Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[br_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd

  Sigma_22_inv <- solve(Sigma_22)
  cross_term <- Sigma_12 %*% Sigma_22_inv %*% t(Sigma_12)
  Sigma_cond <- Sigma_11 - cross_term

  total_dim <- 3 * n_tp + 2
  Sigma <- matrix(0, total_dim, total_dim)
  Sigma[1:(3*n_tp), 1:(3*n_tp)] <- Sigma_11
  Sigma[(3*n_tp+1):total_dim, (3*n_tp+1):total_dim] <- Sigma_22
  Sigma[1:(3*n_tp), (3*n_tp+1):total_dim] <- Sigma_12
  Sigma[(3*n_tp+1):total_dim, 1:(3*n_tp)] <- t(Sigma_12)

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

# ============================================================================
# PARTICIPANT GENERATION
# ============================================================================

generate_participant_twostage <- function(sigma_parts, idx) {
  Sigma_22 <- sigma_parts$Sigma_22
  Sigma_cond <- sigma_parts$Sigma_cond

  stage1 <- mvrnorm(1, mu = c(biomarker_mean, baseline_mean), Sigma = Sigma_22)
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

# ============================================================================
# DESIGN STRUCTURES
# ============================================================================

measurement_weeks_hybrid <- c(4, 8, 9, 10, 11, 12, 16, 20)
measurement_weeks_ol_bdc <- c(4, 8, 12, 16, 17, 18, 19, 20)
measurement_weeks_crossover <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)
measurement_weeks_ol <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)

create_ol_bdc_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        week <= 16 ~ 1,
        week == 17 ~ 1,
        week == 18 & path %in% c(1, 2) ~ 1,
        week == 18 & path %in% c(3, 4) ~ 0,
        week %in% c(19, 20) ~ 0,
        TRUE ~ 1
      ),
      expectancy = if_else(week <= 16, 1.0, 0.5)
    )
}

create_hybrid_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
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

create_crossover_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:2, length.out = n_participants))
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        week <= 10 & path == 1 ~ 1,
        week <= 10 & path == 2 ~ 0,
        week > 10 & path == 1 ~ 0,
        week > 10 & path == 2 ~ 1,
        TRUE ~ NA_real_
      ),
      expectancy = 0.5
    )
}

create_ol_design <- function(n_participants, measurement_weeks) {
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = 1,
      treatment = 1,
      expectancy = 1.0
    )
}

# ============================================================================
# BUILD SIGMA CACHE
# ============================================================================

cat("Building sigma cache...\n")

sigma_structures <- list(

list(design = "hybrid", weeks = measurement_weeks_hybrid, c_bm = 0.3),
  list(design = "hybrid", weeks = measurement_weeks_hybrid, c_bm = 0),
  list(design = "ol_bdc", weeks = measurement_weeks_ol_bdc, c_bm = 0.3),
  list(design = "ol_bdc", weeks = measurement_weeks_ol_bdc, c_bm = 0),
  list(design = "crossover", weeks = measurement_weeks_crossover, c_bm = 0.3),
  list(design = "crossover", weeks = measurement_weeks_crossover, c_bm = 0),
  list(design = "ol", weeks = measurement_weeks_ol, c_bm = 0.3),
  list(design = "ol", weeks = measurement_weeks_ol, c_bm = 0)
)

sigma_cache <- list()
for (i in seq_along(sigma_structures)) {
  sigma_def <- sigma_structures[[i]]
  cache_key <- paste(sigma_def$design, sigma_def$c_bm, sep = "_")
  sigma_cache[[cache_key]] <- list(
    sigma_parts = build_sigma_guaranteed_pd(
      weeks = sigma_def$weeks,
      c.bm = sigma_def$c_bm,
      params = list()
    ),
    weeks = sigma_def$weeks
  )
}

get_cached_sigma <- function(design, c_bm, cache) {
  cache_key <- paste(design, c_bm, sep = "_")
  cache[[cache_key]]
}

cat("Sigma cache ready\n\n")

# ============================================================================
# SINGLE ITERATION FUNCTION (with censoring)
# ============================================================================

run_single_iteration <- function(iter, params, sigma_cache, model_engine) {
  set.seed(iter * 1000 + params$condition_id)

  cached <- get_cached_sigma(params$design, params$biomarker_correlation, sigma_cache)
  sigma_parts <- cached$sigma_parts
  measurement_weeks <- cached$weeks
  idx <- sigma_parts$indices
  effective_bm_corr <- sigma_parts$effective_c.bm

  local_n <- params$n_participants

  # Create design
  if (params$design == "hybrid") {
    trial_design <- create_hybrid_design(local_n, measurement_weeks)
  } else if (params$design == "crossover") {
    trial_design <- create_crossover_design(local_n, measurement_weeks)
  } else if (params$design == "ol") {
    trial_design <- create_ol_design(local_n, measurement_weeks)
  } else {
    trial_design <- create_ol_bdc_design(local_n, measurement_weeks)
  }

  # Generate participant data
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

  bm_mod <- params$biomarker_moderation
  carryover_t_half <- params$carryover_halflife
  beta_0 <- params$beta_0
  beta_1 <- params$beta_1

  # Build trial data
  trial_data <- trial_design |>
    group_by(participant_id) |>
    mutate(timepoint_idx = row_number()) |>
    ungroup() |>
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
      bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
      effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered),
      br_on_drug = weeks_on_drug * effective_BR_rate,
      br_at_discont = if_else(just_discontinued, lag(br_on_drug, default = 0), NA_real_),
      br_at_discont = zoo::na.locf(br_at_discont, na.rm = FALSE),
      carryover_effect = case_when(
        treatment == 1 ~ 1,
        treatment == 0 & carryover_t_half == 0 ~ 0,
        treatment == 0 & carryover_t_half > 0 ~ (1/2)^(tsd / carryover_t_half),
        TRUE ~ 0
      ),
      br_off_drug = if_else(
        treatment == 0 & !is.na(br_at_discont),
        br_at_discont * carryover_effect,
        0
      ),
      BR_mean = if_else(treatment == 1, br_on_drug, br_off_drug),
      weeks_at_discont = if_else(just_discontinued, lag(weeks_on_drug, default = 0), NA_real_),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE),
      effective_drug_weeks = case_when(
        treatment == 1 ~ weeks_on_drug,
        treatment == 0 & !is.na(weeks_at_discont) ~ weeks_at_discont * carryover_effect,
        TRUE ~ 0
      ),
      BR = BR_mean + br_random,
      ER_mean = weeks_with_expectancy * ER_rate,
      ER = ER_mean + er_random,
      TR_mean = weeks_in_trial * TR_rate,
      TR = TR_mean + tr_random,
      response = baseline + BR + ER + TR
    ) |>
    ungroup()

  # Apply censoring (Hendrickson-style: dropout probability per timepoint)
  if (beta_0 > 0 || beta_1 > 0) {
    # Add tracking columns
    trial_data$dropped_out <- FALSE

    # Process each participant sequentially to simulate dropout cascade
    for (pid in unique(trial_data$participant_id)) {
      pid_rows <- which(trial_data$participant_id == pid)
      baseline_val <- trial_data$baseline[pid_rows[1]]
      dropped <- FALSE

      for (row_idx in pid_rows) {
        if (dropped) {
          trial_data$dropped_out[row_idx] <- TRUE
        } else {
          response_val <- trial_data$response[row_idx]
          if (apply_censoring(response_val, baseline_val, beta_0, beta_1)) {
            dropped <- TRUE
            trial_data$dropped_out[row_idx] <- TRUE
          }
        }
      }
    }

    # Remove dropped out observations
    trial_data <- trial_data[!trial_data$dropped_out, ]
    trial_data$dropped_out <- NULL
  }

  # Standardize biomarker
  trial_data <- trial_data |>
    mutate(bm_centered = (biomarker - mean(biomarker)) / sd(biomarker))

  # Count remaining observations
  n_obs <- nrow(trial_data)
  n_participants_remaining <- length(unique(trial_data$participant_id))

  # Fit model
  model_result <- tryCatch({
    if (params$design == "ol") {
      model <- lme(
        response ~ week * bm_centered,
        random = ~ 1 | participant_id,
        correlation = corCAR1(form = ~ week | participant_id),
        data = trial_data,
        control = lmeControl(opt = "optim")
      )
      interaction_term <- "week:bm_centered"
    } else {
      model <- lme(
        response ~ effective_drug_weeks * bm_centered + week,
        random = ~ 1 | participant_id,
        correlation = corCAR1(form = ~ week | participant_id),
        data = trial_data,
        control = lmeControl(opt = "optim")
      )
      interaction_term <- "effective_drug_weeks:bm_centered"
    }

    coefs <- summary(model)$tTable
    idx_coef <- which(rownames(coefs) == interaction_term)

    tibble(
      iteration = iter,
      design = params$design,
      n_participants = local_n,
      biomarker_moderation = params$biomarker_moderation,
      biomarker_correlation = effective_bm_corr,
      carryover_halflife = params$carryover_halflife,
      censoring_type = params$censoring_type,
      beta_0 = beta_0,
      beta_1 = beta_1,
      n_obs = n_obs,
      n_participants_remaining = n_participants_remaining,
      effect_size = coefs[idx_coef, "Value"],
      se = coefs[idx_coef, "Std.Error"],
      t_value = coefs[idx_coef, "t-value"],
      p_value = coefs[idx_coef, "p-value"],
      significant = coefs[idx_coef, "p-value"] < 0.05
    )
  }, error = function(e) {
    tibble(
      iteration = iter,
      design = params$design,
      n_participants = local_n,
      biomarker_moderation = params$biomarker_moderation,
      biomarker_correlation = effective_bm_corr,
      carryover_halflife = params$carryover_halflife,
      censoring_type = params$censoring_type,
      beta_0 = beta_0,
      beta_1 = beta_1,
      n_obs = n_obs,
      n_participants_remaining = n_participants_remaining,
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
# RUN SIMULATION
# ============================================================================

log_file <- "../output/simulation_censoring_log.txt"
if (!dir.exists("../output")) dir.create("../output", recursive = TRUE)
sink(log_file, split = TRUE)

cat("Running CENSORING SENSITIVITY ANALYSIS...\n")
cat("Iterations per condition =", n_iterations, "\n")
cat("Total conditions =", nrow(param_grid), "\n\n")

if (use_parallel && n_cores > 1) {
  plan(multisession, workers = n_cores)
  cat(sprintf("Parallel backend: %d workers\n\n", n_cores))
} else {
  plan(sequential)
  cat("Sequential mode\n\n")
}

start_time <- Sys.time()
results_list <- list()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf(
    "Condition %d/%d: %s, N=%d, bm=%.1f, t½=%.1f, cens=%s\n",
    i, nrow(param_grid), params$design, params$n_participants,
    params$biomarker_moderation, params$carryover_halflife, params$censoring_type
  ))

  local_params <- params
  local_sigma_cache <- sigma_cache
  local_model_engine <- model_engine

  condition_results <- future_map_dfr(
    1:n_iterations,
    function(iter) {
      run_single_iteration(iter, local_params, local_sigma_cache, local_model_engine)
    },
    .options = furrr_options(
      seed = TRUE,
      globals = list(
        run_single_iteration = run_single_iteration,
        get_cached_sigma = get_cached_sigma,
        generate_participant_twostage = generate_participant_twostage,
        create_hybrid_design = create_hybrid_design,
        create_ol_bdc_design = create_ol_bdc_design,
        create_crossover_design = create_crossover_design,
        create_ol_design = create_ol_design,
        apply_censoring = apply_censoring,
        local_params = local_params,
        local_sigma_cache = local_sigma_cache,
        local_model_engine = local_model_engine,
        biomarker_mean = biomarker_mean,
        baseline_mean = baseline_mean,
        biomarker_sd = biomarker_sd,
        BR_rate = BR_rate,
        ER_rate = ER_rate,
        TR_rate = TR_rate
      ),
      packages = c("tidyverse", "nlme", "lme4", "MASS", "zoo")
    )
  )

  results_list[[i]] <- condition_results

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate <- i / elapsed
  remaining <- (nrow(param_grid) - i) / rate
  cat(sprintf("  %.1f sec elapsed, ETA: %.0f sec\n", elapsed, remaining))
}

results <- bind_rows(results_list)
plan(sequential)

# ============================================================================
# SUMMARIZE RESULTS
# ============================================================================

summary_results <- results |>
  group_by(design, n_participants, biomarker_moderation, carryover_halflife,
           censoring_type) |>
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    mean_n_obs = mean(n_obs, na.rm = TRUE),
    mean_n_remaining = mean(n_participants_remaining, na.rm = TRUE),
    pct_retained = mean(n_participants_remaining / n_participants * 100, na.rm = TRUE),
    n = n(),
    n_errors = sum(is.na(significant)),
    .groups = "drop"
  ) |>
  mutate(
    n_success = round(power * n),
    ci_lower = (power + 1.96^2 / (2 * n) -
      1.96 * sqrt((power * (1 - power) + 1.96^2 / (4 * n)) / n)) /
      (1 + 1.96^2 / n),
    ci_upper = (power + 1.96^2 / (2 * n) +
      1.96 * sqrt((power * (1 - power) + 1.96^2 / (4 * n)) / n)) /
      (1 + 1.96^2 / n),
    ci_lower = pmax(0, ci_lower),
    ci_upper = pmin(1, ci_upper)
  )

print(summary_results, n = 50)

# ============================================================================
# VISUALIZATION - Hendrickson Figure 4 Style
# ============================================================================

library(viridis)

plot_data <- summary_results |>
  filter(!is.na(power)) |>
  mutate(
    censoring_label = factor(
      censoring_type,
      levels = c("high_dropout", "more_biased", "more_flat", "balanced", "none"),
      labels = c("High dropout", "More biased", "More flat", "Balanced", "No censoring")
    ),
    bm_label = factor(
      biomarker_moderation,
      levels = c(0, 0.3, 0.6)
    ),
    carryover_label = factor(
      paste0("t½=", carryover_halflife),
      levels = paste0("t½=", c(0, 0.1, 0.2))
    ),
    design_label = factor(
      case_when(
        design == "crossover" ~ "Crossover",
        design == "hybrid" ~ "N-of-1",
        design == "ol_bdc" ~ "OL+BDC",
        design == "ol" ~ "Open-Label"
      ),
      levels = c("Open-Label", "OL+BDC", "Crossover", "N-of-1")
    ),
    n_label = factor(paste0("N=", n_participants), levels = c("N=35", "N=70")),
    power_pct = round(power * 100)
  )

# Hendrickson Figure 4A style: Effect of biomarker and censoring on power
# Layout: Carryover on x-axis, censoring on y-axis, faceted by design and N

p_main <- ggplot(
  plot_data,
  aes(x = carryover_label, y = censoring_label, fill = power)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", power)),
    color = "black", size = 2.2
  ) +
  scale_fill_gradient2(
    name = "Power",
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1)
  ) +
  facet_grid(
    n_label + bm_label ~ design_label,
    labeller = labeller(
      n_label = label_value,
      bm_label = function(x) paste0("BM=", x),
      design_label = label_value
    )
  ) +
  labs(
    title = "Effect of Trial Design, Censoring, and Carryover on Power",
    subtitle = sprintf(
      "%d iterations | Hendrickson-style censoring sensitivity analysis",
      n_iterations
    ),
    x = "Carryover Half-Life",
    y = "Censoring Parameters"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 8),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  )

print(p_main)

# Save outputs
output_dir <- "../output"
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

ggsave(
  file.path(output_dir, sprintf("power_censoring_%s.pdf", timestamp)),
  p_main,
  width = 14, height = 16
)
ggsave(
  file.path(output_dir, sprintf("power_censoring_%s.png", timestamp)),
  p_main,
  width = 14, height = 16, dpi = 300
)

save(
  results, summary_results, param_grid, censoring_params, n_iterations,
  file = file.path(output_dir, "simulation_censoring_results.RData")
)

cat("\n", strrep("=", 70), "\n")
cat("CENSORING SENSITIVITY ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("Results saved to:", output_dir, "\n")
cat(sprintf("  - power_censoring_%s.pdf/png\n", timestamp))
cat("  - simulation_censoring_results.RData\n")

sink()
