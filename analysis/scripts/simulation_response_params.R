# ============================================================================
# RESPONSE PARAMETERS SENSITIVITY ANALYSIS
# ============================================================================
# Analogous to Hendrickson Figure 5: Power as a function of response trajectory
# parameters.
#
# Our model uses a LINEAR response model (not Gompertz):
#   BR = weeks_on_drug * BR_rate * (1 + biomarker_moderation * bm_centered)
#   ER = week * ER_rate * expectancy
#   TR = week * TR_rate
#
# This script varies the slope parameters (BR_rate, ER_rate, TR_rate) to
# examine how power depends on the relative magnitudes of these effects.
#
# Usage: Rscript simulation_response_params.R
#        EXPLORATORY_MODE=TRUE Rscript simulation_response_params.R
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

n_cores <- max(1, parallel::detectCores() - 1)
use_parallel <- TRUE
model_engine <- "lme"

# Base response parameters (will be varied in sensitivity analysis)
BR_rate_base <- 0.5
ER_rate_base <- 0.2
TR_rate_base <- 0.1

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

carryover_mode <- "halflife"

# ============================================================================
# RESPONSE PARAMETER GRID (Hendrickson Figure 5B style)
# ============================================================================
# Full factorial design: BR_rate × ER_rate × TR_rate
# Each varied at 3 levels to produce a range of power results
#
# Heatmap layout:
#   - Y-axis: BR_rate
#   - X-axis: ER_rate
#   - Facet rows: TR_rate
#   - Facet columns: Trial design

rate_levels <- c(0.05, 0.15, 0.3)

response_params_grid <- expand_grid(
  BR_rate = rate_levels,
  ER_rate = rate_levels,
  TR_rate = rate_levels
) |>
  mutate(
    param_set = sprintf("BR%.2f_ER%.2f_TR%.2f", BR_rate, ER_rate, TR_rate)
  )

# ============================================================================
# PARAMETER GRID - Response Parameters Sensitivity
# ============================================================================
# Full factorial grid:
# - 4 designs
# - 2 sample sizes
# - 27 response parameter sets (3×3×3)
# - Fixed: biomarker_moderation=0.3, carryover=0, no censoring
# Total: 4 × 2 × 27 = 216 conditions

param_grid <- expand_grid(
  design = c("ol", "ol_bdc", "crossover", "hybrid"),
  n_participants = c(35, 70),
  param_set = response_params_grid$param_set
) |>
  left_join(response_params_grid, by = "param_set") |>
  mutate(
    biomarker_moderation = 0.3,
    biomarker_correlation = 0.3,
    carryover_halflife = 0,
    condition_id = row_number()
  )

cat("Total conditions:", nrow(param_grid), "\n")
cat("Conditions breakdown:\n")
cat("  Designs:", length(unique(param_grid$design)), "\n")
cat("  Sample sizes:", paste(unique(param_grid$n_participants), collapse = ", "), "\n")
cat("  Response parameter sets:", length(unique(param_grid$param_set)), "\n\n")

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

  Sigma_12 <- matrix(0, 3 * n_tp, 2)
  Sigma_12[br_idx, 1] <- effective_c.bm * within_subject_sd * biomarker_sd
  Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[br_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd

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
# DESIGN STRUCTURES
# ============================================================================

measurement_weeks_hybrid <- c(4, 8, 9, 10, 11, 12, 16, 20)
measurement_weeks_ol_bdc <- c(4, 8, 12, 16, 17, 18, 19, 20)
measurement_weeks_crossover <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)
measurement_weeks_ol <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)

create_ol_bdc_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))
  map_dfr(1:n_participants, function(pid) {
    path <- path_assignment[pid]
    tibble(
      participant_id = pid,
      week = measurement_weeks,
      timepoint = 1:length(measurement_weeks),
      on_drug = c(rep(TRUE, 4), path %in% c(1, 2), FALSE, FALSE, FALSE),
      blinded = c(rep(FALSE, 4), TRUE, TRUE, TRUE, TRUE),
      expectancy = c(rep(1, 4), rep(0.5, 4)),
      randomization_path = path
    )
  })
}

create_hybrid_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))
  map_dfr(1:n_participants, function(pid) {
    path <- path_assignment[pid]
    bdc_on_drug <- path %in% c(1, 2)
    cross1_drug <- path %in% c(1, 3)
    cross2_drug <- !cross1_drug
    tibble(
      participant_id = pid,
      week = measurement_weeks,
      timepoint = 1:length(measurement_weeks),
      on_drug = c(TRUE, TRUE, bdc_on_drug, FALSE, cross1_drug, cross1_drug,
                  cross2_drug, cross2_drug),
      blinded = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      expectancy = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
      randomization_path = path
    )
  })
}

create_crossover_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:2, length.out = n_participants))
  map_dfr(1:n_participants, function(pid) {
    path <- path_assignment[pid]
    drug_first <- path == 1
    tibble(
      participant_id = pid,
      week = measurement_weeks,
      timepoint = 1:length(measurement_weeks),
      on_drug = if (drug_first) c(rep(TRUE, 4), rep(FALSE, 4))
                else c(rep(FALSE, 4), rep(TRUE, 4)),
      blinded = TRUE,
      expectancy = 0.5,
      randomization_path = path
    )
  })
}

create_ol_design <- function(n_participants, measurement_weeks) {
  map_dfr(1:n_participants, function(pid) {
    tibble(
      participant_id = pid,
      week = measurement_weeks,
      timepoint = 1:length(measurement_weeks),
      on_drug = TRUE,
      blinded = FALSE,
      expectancy = 1,
      randomization_path = 1
    )
  })
}

# ============================================================================
# RESPONSE GENERATION (with variable rates)
# ============================================================================

generate_response <- function(design_df, params, sigma_parts, idx,
                              BR_rate, ER_rate, TR_rate) {
  n_participants <- max(design_df$participant_id)
  carryover_halflife <- params$carryover_halflife
  biomarker_moderation <- params$biomarker_moderation

  map_dfr(1:n_participants, function(pid) {
    p_data <- generate_participant_twostage(sigma_parts, idx)
    p_design <- design_df |> filter(participant_id == pid)

    last_drug_week <- NA
    responses <- numeric(nrow(p_design))

    for (i in seq_len(nrow(p_design))) {
      row <- p_design[i, ]
      week <- row$week
      tp <- row$timepoint

      br_component <- 0
      er_component <- 0
      tr_component <- 0

      tr_component <- TR_rate * week + p_data$tr_random[tp]

      if (row$on_drug) {
        br_base <- BR_rate * week
        bm_centered <- (p_data$biomarker - biomarker_mean) / biomarker_sd
        br_modulation <- 1 + biomarker_moderation * bm_centered
        br_component <- br_base * br_modulation + p_data$br_random[tp]
        last_drug_week <- week
      } else if (!is.na(last_drug_week) && carryover_halflife > 0) {
        time_since_drug <- week - last_drug_week
        carryover_prop <- calculate_carryover(time_since_drug, carryover_halflife)
        br_base <- BR_rate * last_drug_week
        bm_centered <- (p_data$biomarker - biomarker_mean) / biomarker_sd
        br_modulation <- 1 + biomarker_moderation * bm_centered
        br_component <- br_base * br_modulation * carryover_prop +
                       p_data$br_random[tp] * carryover_prop
      }

      er_component <- ER_rate * week * row$expectancy + p_data$er_random[tp]

      response <- p_data$baseline - (br_component + er_component + tr_component)
      responses[i] <- response
    }

    p_design |>
      mutate(
        biomarker = p_data$biomarker,
        baseline = p_data$baseline,
        response = responses,
        censored = FALSE
      )
  })
}

# ============================================================================
# MODEL FITTING
# ============================================================================

fit_model_safe <- function(trial_data, design_type) {
  tryCatch({
    if (design_type == "ol") {
      model <- lme(
        response ~ week + biomarker + week:biomarker,
        random = ~ 1 | participant_id,
        data = trial_data,
        correlation = corCAR1(form = ~ week | participant_id),
        control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      coef_name <- "week:biomarker"
    } else {
      model <- lme(
        response ~ on_drug + biomarker + on_drug:biomarker,
        random = ~ 1 | participant_id,
        data = trial_data,
        correlation = corCAR1(form = ~ week | participant_id),
        control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      )
      coef_name <- "on_drugTRUE:biomarker"
    }

    coef_table <- summary(model)$tTable
    if (coef_name %in% rownames(coef_table)) {
      p_value <- coef_table[coef_name, "p-value"]
      estimate <- coef_table[coef_name, "Value"]
      std_error <- coef_table[coef_name, "Std.Error"]
    } else {
      p_value <- NA
      estimate <- NA
      std_error <- NA
    }

    list(
      p_value = p_value,
      estimate = estimate,
      std_error = std_error,
      converged = TRUE,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      p_value = NA,
      estimate = NA,
      std_error = NA,
      converged = FALSE,
      error = conditionMessage(e)
    )
  })
}

# ============================================================================
# SINGLE ITERATION RUNNER
# ============================================================================

run_single_iteration <- function(params, sigma_cache) {
  design_type <- params$design
  n_participants <- params$n_participants
  BR_rate <- params$BR_rate
  ER_rate <- params$ER_rate
  TR_rate <- params$TR_rate

  measurement_weeks <- switch(
    design_type,
    "hybrid" = measurement_weeks_hybrid,
    "ol_bdc" = measurement_weeks_ol_bdc,
    "crossover" = measurement_weeks_crossover,
    "ol" = measurement_weeks_ol
  )

  cache_key <- paste(design_type, params$biomarker_moderation, sep = "_")
  sigma_parts <- sigma_cache[[cache_key]]
  idx <- sigma_parts$indices

  design_df <- switch(
    design_type,
    "hybrid" = create_hybrid_design(n_participants, measurement_weeks),
    "ol_bdc" = create_ol_bdc_design(n_participants, measurement_weeks),
    "crossover" = create_crossover_design(n_participants, measurement_weeks),
    "ol" = create_ol_design(n_participants, measurement_weeks)
  )

  trial_data <- generate_response(design_df, params, sigma_parts, idx,
                                  BR_rate, ER_rate, TR_rate)

  fit_result <- fit_model_safe(trial_data, design_type)

  tibble(
    design = design_type,
    n_participants = n_participants,
    param_set = params$param_set,
    BR_rate = BR_rate,
    ER_rate = ER_rate,
    TR_rate = TR_rate,
    biomarker_moderation = params$biomarker_moderation,
    carryover_halflife = params$carryover_halflife,
    p_value = fit_result$p_value,
    estimate = fit_result$estimate,
    std_error = fit_result$std_error,
    converged = fit_result$converged
  )
}

# ============================================================================
# BUILD SIGMA CACHE
# ============================================================================

cat("Building sigma cache...\n")

unique_conditions <- param_grid |>
  select(design, biomarker_moderation, biomarker_correlation) |>
  distinct()

sigma_cache <- list()

for (i in seq_len(nrow(unique_conditions))) {
  cond <- unique_conditions[i, ]
  design_type <- cond$design
  bm_mod <- cond$biomarker_moderation
  bm_corr <- cond$biomarker_correlation

  measurement_weeks <- switch(
    design_type,
    "hybrid" = measurement_weeks_hybrid,
    "ol_bdc" = measurement_weeks_ol_bdc,
    "crossover" = measurement_weeks_crossover,
    "ol" = measurement_weeks_ol
  )

  cache_key <- paste(design_type, bm_mod, sep = "_")
  sigma_cache[[cache_key]] <- build_sigma_guaranteed_pd(
    measurement_weeks, bm_corr, list(biomarker_moderation = bm_mod)
  )
}

cat("Sigma cache ready\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

cat("Running RESPONSE PARAMETERS SENSITIVITY ANALYSIS...\n")
cat("Iterations per condition =", n_iterations, "\n")
cat("Total conditions =", nrow(param_grid), "\n\n")

if (use_parallel && n_cores > 1) {
  plan(multisession, workers = n_cores)
  cat("Parallel backend:", n_cores, "workers\n\n")
} else {
  plan(sequential)
  cat("Running sequentially\n\n")
}

results <- list()
start_time <- Sys.time()

for (cond_idx in seq_len(nrow(param_grid))) {
  params <- param_grid[cond_idx, ]
  iter_start <- Sys.time()

  iter_results <- future_map_dfr(
    1:n_iterations,
    function(iter) {
      run_single_iteration(params, sigma_cache)
    },
    .options = furrr_options(seed = TRUE)
  )

  results[[cond_idx]] <- iter_results |>
    mutate(iteration = 1:n_iterations)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  eta <- elapsed / cond_idx * (nrow(param_grid) - cond_idx)

  cat(sprintf(
    "Condition %d/%d: %s, N=%d, %s\n",
    cond_idx, nrow(param_grid),
    params$design, params$n_participants, params$param_set
  ))
  cat(sprintf("  %.1f sec elapsed, ETA: %.0f sec\n", elapsed, eta))
}

results <- bind_rows(results)

# ============================================================================
# SUMMARIZE RESULTS
# ============================================================================

summary_results <- results |>
  group_by(design, n_participants, param_set, BR_rate, ER_rate, TR_rate,
           biomarker_moderation, carryover_halflife) |>
  summarise(
    power = mean(p_value < 0.05, na.rm = TRUE),
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE),
    convergence_rate = mean(converged, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop"
  )

# ============================================================================
# SAVE RESULTS
# ============================================================================

output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save(
  results, summary_results, param_grid, response_params_grid, n_iterations,
  file = file.path(output_dir, "simulation_response_params_results.RData")
)

cat("\nResults saved to:", file.path(output_dir,
    "simulation_response_params_results.RData"), "\n")

# ============================================================================
# VISUALIZATION - Hendrickson Figure 5B Style
# ============================================================================
# Layout: BR_rate on Y-axis, ER_rate on X-axis
# Facet grid: TR_rate (rows) × design (columns)

plot_data <- summary_results |>
  filter(!is.na(power)) |>
  mutate(
    design_label = factor(
      case_when(
        design == "crossover" ~ "CO",
        design == "hybrid" ~ "N-of-1",
        design == "ol_bdc" ~ "OL+BDC",
        design == "ol" ~ "OL"
      ),
      levels = c("OL", "OL+BDC", "CO", "N-of-1")
    ),
    BR_label = factor(BR_rate, levels = rev(rate_levels)),
    ER_label = factor(ER_rate, levels = rate_levels),
    TR_label = factor(paste0("TR=", TR_rate), levels = paste0("TR=", rate_levels)),
    power_display = sprintf("%.2f", power)
  )

# Function to create plot for a single sample size
make_rate_heatmap <- function(data, n_val) {
  plot_subset <- data |> filter(n_participants == n_val)

  ggplot(
    plot_subset,
    aes(x = ER_label, y = BR_label, fill = power)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = power_display), color = "black", size = 2.2) +
    scale_fill_gradient2(
      name = "Power",
      low = "#d73027", mid = "#fee08b", high = "#1a9850",
      midpoint = 0.5, limits = c(0, 1)
    ) +
    facet_grid(TR_label ~ design_label) +
    labs(
      title = sprintf("N = %d", n_val),
      x = "ER rate",
      y = "BR rate"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "right",
      panel.spacing = unit(0.3, "lines")
    )
}

p_n35 <- make_rate_heatmap(plot_data, 35)
p_n70 <- make_rate_heatmap(plot_data, 70)

p_heatmap <- p_n35 / p_n70 +
  plot_annotation(
    title = "Effect of Response Parameter Rate Values on Power",
    subtitle = sprintf(
      "%d iterations | BR/ER/TR rates: %s | biomarker_mod=0.3",
      n_iterations, paste(rate_levels, collapse = ", ")
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  )

# ============================================================================
# SAVE PLOTS
# ============================================================================

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

ggsave(
  file.path(output_dir, sprintf("power_response_params_%s.pdf", timestamp)),
  p_heatmap,
  width = 10, height = 8
)
ggsave(
  file.path(output_dir, sprintf("power_response_params_%s.png", timestamp)),
  p_heatmap,
  width = 10, height = 8, dpi = 300
)

cat("Plots saved to:", output_dir, "\n")
cat(sprintf("  - power_response_params_%s.pdf/png\n", timestamp))

# Print summary table (condensed)
cat("\n=== POWER SUMMARY (mean across TR levels) ===\n")
summary_results |>
  group_by(design, n_participants, BR_rate, ER_rate) |>
  summarise(mean_power = mean(power, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = design, values_from = mean_power) |>
  arrange(n_participants, BR_rate, ER_rate) |>
  print(n = 30)
