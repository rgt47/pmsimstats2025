# N-of-1 Trial Simulation - Evenly-Spaced Measurement Designs
# OL, Crossover, and Parallel designs with 4 evenly-spaced measurement points
# These designs have uniform temporal sampling

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(MASS)
  library(conflicted)
})

# Resolve conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(nlme::lme)

# ============================================================================
# PARAMETERS
# ============================================================================

n_participants <- 70
n_iterations <- 20

# Three-factor response model - RATE-BASED (points per week)
BR_rate <- 0.5  # Biological Response: drug improvement rate
ER_rate <- 0.2  # Expectancy Response: placebo improvement rate
TR_rate <- 0.1  # Time-variant Response: natural improvement rate
treatment_effect <- BR_rate  # Alias for display purposes

baseline_mean <- 10.0       # Mean baseline response
between_subject_sd <- 2.0   # SD of participant random effects
within_subject_sd <- 1.8    # SD of measurement noise
biomarker_mean <- 5.0       # Mean biomarker value
biomarker_sd <- 2.0         # SD of biomarker

# Correlation parameters (from Hendrickson)
c.br <- 0.8  # BR autocorrelation across timepoints
c.er <- 0.8  # ER autocorrelation across timepoints
c.tr <- 0.8  # TR autocorrelation across timepoints
c.cf1t <- 0.2  # Same-time cross-correlation
c.cfct <- 0.1  # Different-time cross-correlation
c.bm_baseline <- 0.3   # Biomarker-baseline correlation
c.baseline_resp <- 0.4 # Baseline-response correlation

# Allowed correlation values (finite grid for consistent results)
allowed_correlations <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

# Parameter grid - OL, Crossover, and Parallel designs only
param_grid <- bind_rows(
  # Open-label design (power conditions)
  expand_grid(
    design = "ol",
    biomarker_moderation = c(0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    carryover = c(0)
  ),
  # Type I error condition for OL
  expand_grid(
    design = "ol",
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    carryover = c(0)
  ),
  # Crossover design (power conditions)
  expand_grid(
    design = "crossover",
    biomarker_moderation = c(0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    carryover = c(0, 0.5, 1)
  ),
  # Type I error condition for Crossover
  expand_grid(
    design = "crossover",
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    carryover = c(0, 0.5, 1)
  ),
  # Parallel design (power conditions)
  expand_grid(
    design = "parallel",
    biomarker_moderation = c(0, 0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    carryover = c(0)
  )
)

# ============================================================================
# PRE-SIMULATION PARAMETER VALIDATION
# ============================================================================

cat("\n")
cat("=" %+% strrep("=", 79) %+% "\n")
cat("PRE-SIMULATION PARAMETER VALIDATION\n")
cat("=" %+% strrep("=", 79) %+% "\n\n")

# Create response and baseline parameters for validation
resp_param <- tibble(
  cat = c("biological_response", "expectancy_response", "time_variant_response"),
  max = c(1.0, 1.0, 1.0),
  sd = c(within_subject_sd, within_subject_sd, within_subject_sd)
)

baseline_param <- tibble(
  cat = c("biomarker", "baseline"),
  m = c(biomarker_mean, baseline_mean),
  sd = c(biomarker_sd, between_subject_sd)
)

# Create model parameters for validation
model_params <- list(
  N = n_participants,
  c.br = c.br,
  c.er = c.er,
  c.tr = c.tr,
  c.cf1t = c.cf1t,
  c.cfct = c.cfct,
  c.bm_baseline = c.bm_baseline,
  c.baseline_resp = c.baseline_resp
)

# Validate for each design separately (different designs may have different behavior)
# All evenly-spaced designs use the same measurement schedule
designs_to_validate <- list(
  ol = measurement_weeks,
  crossover = measurement_weeks,
  parallel = measurement_weeks
)

all_valid_rows <- tibble()

for (design_name in names(designs_to_validate)) {
  cat(sprintf("Validating %s design (measurement weeks: %s)...\n",
              toupper(design_name),
              paste(designs_to_validate[[design_name]], collapse = ", ")))

  # Create representative trial design for this design
  meas_weeks <- designs_to_validate[[design_name]]

  if (design_name == "ol") {
    trial_design <- create_ol_design(n_participants, meas_weeks)
  } else if (design_name == "crossover") {
    trial_design <- create_crossover_design(n_participants, meas_weeks)
  } else {
    trial_design <- create_parallel_design(n_participants, meas_weeks)
  }

  # Filter param_grid to this design
  design_params <- param_grid %>% filter(design == design_name)

  # Validate parameter grid for this design
  validation_result <- validate_parameter_grid(
    param_grid = design_params,
    trial_design = trial_design,
    model_params = model_params,
    resp_param = resp_param,
    baseline_param = baseline_param,
    verbose = FALSE
  )

  # Report validation results
  cat(sprintf("  Valid combinations:   %d / %d (%.1f%%)\n",
              validation_result$n_valid,
              nrow(design_params),
              100 * validation_result$n_valid / nrow(design_params)))

  if (validation_result$n_invalid > 0) {
    cat(sprintf("  Invalid combinations: %d\n\n", validation_result$n_invalid))
    cat("  Details of invalid combinations:\n")
    print(validation_result$invalid_combinations)
    cat("\n  Reasons for failure:\n")
    for (reason in validation_result$invalid_reasons) {
      cat(sprintf("    • %s\n", reason))
    }
    cat("\n")
  }

  # Report condition numbers
  if (length(validation_result$condition_numbers) > 0) {
    cond <- validation_result$condition_numbers
    cat(sprintf("  Condition number (κ) statistics:\n"))
    cat(sprintf("    Mean:   %.1f\n", mean(cond)))
    cat(sprintf("    Median: %.1f\n", median(cond)))
    cat(sprintf("    Min:    %.1f (best conditioned)\n", min(cond)))
    cat(sprintf("    Max:    %.1f (worst conditioned)\n", max(cond)))

    if (max(cond) > 100) {
      cat(sprintf("    ⚠ WARNING: Some matrices are ill-conditioned (κ > 100)\n"))
      cat(sprintf("      You may see convergence warnings during simulation.\n"))
    }
  }
  cat("\n")

  # Add valid combinations to overall list with design label preserved
  all_valid_rows <- bind_rows(
    all_valid_rows,
    validation_result$valid_combinations
  )
}

# Update param_grid to only include valid combinations
original_n <- nrow(param_grid)
param_grid <- all_valid_rows

if (original_n > nrow(param_grid)) {
  cat(sprintf("FILTERED: %d -> %d parameter combinations\n", original_n, nrow(param_grid)))
  cat(sprintf("          %d combinations excluded due to validation failures\n\n", original_n - nrow(param_grid)))
} else {
  cat(sprintf("✓ All %d parameter combinations passed validation\n\n", nrow(param_grid)))
}

cat("=" %+% strrep("=", 79) %+% "\n\n")

# ============================================================================
# BUILD SIGMA WITH GUARANTEED PD (Time-Based AR(1))
# ============================================================================

build_sigma_guaranteed_pd <- function(weeks, c.bm, params) {
  n_tp <- length(weeks)

  br_idx <- 1:n_tp
  er_idx <- (n_tp + 1):(2 * n_tp)
  tr_idx <- (2 * n_tp + 1):(3 * n_tp)
  bm_idx <- 3 * n_tp + 1
  bl_idx <- 3 * n_tp + 2

  # STAGE 1: Build Sigma_22 (2x2) - Biomarker & Baseline
  Sigma_22 <- matrix(c(
    biomarker_sd^2,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    between_subject_sd^2
  ), 2, 2)

  # STAGE 2: Build Sigma_11 (3*n_tp x 3*n_tp) - Response components
  build_ar1_time <- function(weeks, rho, sigma) {
    n <- length(weeks)
    Cov <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        time_lag <- abs(weeks[i] - weeks[j])
        Cov[i, j] <- sigma^2 * rho^time_lag
      }
    }
    return(Cov)
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

  # STAGE 3: Build Sigma_12 (no on-the-fly adjustment - parameters are pre-validated)
  effective_c.bm <- c.bm
  effective_c.baseline <- c.baseline_resp

  Sigma_12 <- matrix(0, 3 * n_tp, 2)
  Sigma_12[br_idx, 1] <- effective_c.bm * within_subject_sd * biomarker_sd
  Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[br_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd

  # Verify positive definiteness (should always pass if parameters were validated correctly)
  Sigma_22_inv <- solve(Sigma_22)
  cross_term <- Sigma_12 %*% Sigma_22_inv %*% t(Sigma_12)
  Sigma_cond <- Sigma_11 - cross_term
  min_eig <- min(eigen(Sigma_cond, only.values = TRUE)$values)

  if (min_eig <= 1e-6) {
    stop(sprintf(
      paste(
        "\nUNEXPECTED: Parameter combination failed positive definiteness check!",
        "This should NOT happen - parameters were pre-validated.",
        "c.bm = %.2f",
        "Min eigenvalue = %.2e",
        "",
        "This indicates either:",
        "  1. A bug in the validate_parameter_grid() function, or",
        "  2. Inconsistent covariance matrix construction between validation and simulation",
        "",
        "Please investigate the mismatch between validate_parameter_grid()",
        "and build_sigma_guaranteed_pd().",
        sep = "\n"
      ),
      effective_c.bm, min_eig))
  }

  # Build full sigma
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

  return(list(
    Sigma = Sigma,
    Sigma_22 = Sigma_22,
    Sigma_cond = Sigma_cond,
    Sigma_12 = Sigma_12,
    indices = indices,
    effective_c.bm = effective_c.bm
  ))
}

# ============================================================================
# TWO-STAGE GENERATION
# ============================================================================

generate_participant_twostage <- function(sigma_parts, idx) {
  Sigma_22 <- sigma_parts$Sigma_22
  Sigma_cond <- sigma_parts$Sigma_cond
  Sigma_12 <- sigma_parts$Sigma_12
  n_tp <- idx$n_tp

  # Stage 1: Generate (biomarker, baseline)
  stage1 <- mvrnorm(1, mu = c(biomarker_mean, baseline_mean), Sigma = Sigma_22)
  biomarker <- stage1[1]
  baseline <- stage1[2]

  # Stage 2: Generate responses conditional on stage 1
  z <- c(biomarker - biomarker_mean, baseline - baseline_mean)
  cond_mean <- as.vector(Sigma_12 %*% solve(Sigma_22) %*% z)
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
# DESIGN STRUCTURES - EVENLY-SPACED SCHEDULE
# ============================================================================

# 4 evenly-spaced points (6-week gaps) for OL, Crossover, Parallel
measurement_weeks <- c(2, 8, 14, 20)

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

create_crossover_design <- function(n_participants, measurement_weeks) {
  sequence_assignment <- sample(rep(1:2, length.out = n_participants))
  midpoint <- median(measurement_weeks)

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = sequence_assignment[participant_id],
      treatment = case_when(
        week <= midpoint & path == 1 ~ 1,
        week <= midpoint & path == 2 ~ 0,
        week > midpoint & path == 1 ~ 0,
        week > midpoint & path == 2 ~ 1,
        TRUE ~ NA_real_
      ),
      expectancy = 0.5
    )
}

create_parallel_design <- function(n_participants, measurement_weeks) {
  treatment_assignment <- sample(rep(0:1, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = treatment_assignment[participant_id] + 1,
      treatment = treatment_assignment[participant_id],
      expectancy = 0.5
    )
}

# ============================================================================
# RUN SIMULATION
# ============================================================================

# Set up logging
log_file <- "../output/simulation_evenly_spaced_log.txt"
if (!dir.exists("../output")) dir.create("../output", recursive = TRUE)
sink(log_file, split = TRUE)

cat("Running EVENLY-SPACED designs simulation...\n")
cat("Designs: OL, Crossover, Parallel (4 measurement points)\n")
cat("Schedule: weeks", paste(measurement_weeks, collapse = ", "), "\n")
cat("Treatment effect =", treatment_effect, "\n")
cat("Participants =", n_participants, "\n")
cat("Iterations per condition =", n_iterations, "\n\n")

results <- tibble()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf(
    "Condition %d: design=%s, bm_mod=%.2f, carryover=%.2f\n",
    i, params$design, params$biomarker_moderation, params$carryover
  ))

  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    if (params$design == "ol") {
      trial_design <- create_ol_design(n_participants, measurement_weeks)
    } else if (params$design == "crossover") {
      trial_design <- create_crossover_design(n_participants, measurement_weeks)
    } else {
      trial_design <- create_parallel_design(n_participants, measurement_weeks)
    }

    n_timepoints <- length(measurement_weeks)

    sigma_parts <- build_sigma_guaranteed_pd(
      measurement_weeks,
      params$biomarker_correlation,
      params
    )
    idx <- sigma_parts$indices
    effective_bm_corr <- sigma_parts$effective_c.bm

    all_participant_data <- list()
    for (pid in 1:n_participants) {
      result <- generate_participant_twostage(sigma_parts, idx)
      all_participant_data[[pid]] <- tibble(
        participant_id = pid,
        biomarker = result$biomarker,
        baseline = result$baseline,
        br_random = result$br_random,
        er_random = result$er_random,
        tr_random = result$tr_random
      )
    }

    participant_data <- bind_rows(all_participant_data) %>%
      group_by(participant_id) %>%
      mutate(timepoint_idx = row_number()) %>%
      ungroup()

    bm_mod <- params$biomarker_moderation

    trial_data <- trial_design %>%
      group_by(participant_id) %>%
      mutate(timepoint_idx = row_number()) %>%
      ungroup() %>%
      left_join(participant_data, by = c("participant_id", "timepoint_idx")) %>%
      group_by(participant_id) %>%
      arrange(week) %>%
      mutate(
        weeks_on_drug = cumsum(treatment),
        weeks_with_expectancy = cumsum(expectancy),
        weeks_in_trial = week - min(week),
        time_off = cumsum(treatment == 0),
        carryover_effect = as.numeric(
          treatment == 0 & lag(treatment, default = 1) == 1
        ),
        bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
        BR_mean = {
          effective_BR_rate <- BR_rate * (1 + bm_mod * bm_centered)
          br_accumulated <- lag(weeks_on_drug, default = 0) * effective_BR_rate
          first_off <- treatment == 0 & lag(treatment, default = 1) == 1
          ifelse(treatment == 1,
            weeks_on_drug * effective_BR_rate,
            ifelse(first_off, br_accumulated * params$carryover, 0)
          )
        },
        BR = BR_mean + br_random,
        ER_mean = weeks_with_expectancy * ER_rate,
        ER = ER_mean + er_random,
        TR_mean = weeks_in_trial * TR_rate,
        TR = TR_mean + tr_random,
        response = baseline + BR + ER + TR
      ) %>%
      ungroup()

    model_result <- tryCatch({
      trial_data <- trial_data %>%
        mutate(bm_centered = biomarker - mean(biomarker))

      if (params$design == "ol") {
        # OL: use time x biomarker interaction
        model <- lme(
          response ~ bm_centered * week,
          random = ~ 1 | participant_id,
          correlation = corCAR1(form = ~ week | participant_id),
          data = trial_data,
          control = lmeControl(opt = "optim")
        )
        interaction_term <- "bm_centered:week"
      } else {
        # Crossover/Parallel: use treatment x biomarker interaction
        if (params$carryover > 0) {
          model <- lme(
            response ~ treatment * bm_centered + week + carryover_effect,
            random = ~ 1 | participant_id,
            correlation = corCAR1(form = ~ week | participant_id),
            data = trial_data,
            control = lmeControl(opt = "optim")
          )
        } else {
          model <- lme(
            response ~ treatment * bm_centered + week,
            random = ~ 1 | participant_id,
            correlation = corCAR1(form = ~ week | participant_id),
            data = trial_data,
            control = lmeControl(opt = "optim")
          )
        }
        interaction_term <- "treatment:bm_centered"
      }

      coefs <- summary(model)$tTable
      idx_coef <- which(rownames(coefs) == interaction_term)

      tibble(
        iteration = iter,
        design = params$design,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,
        carryover = params$carryover,
        effect_size = coefs[idx_coef, "Value"],
        se = coefs[idx_coef, "Std.Error"],
        t_value = coefs[idx_coef, "t-value"],
        p_value = coefs[idx_coef, "p-value"],
        significant = coefs[idx_coef, "p-value"] < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(
        iteration = iter,
        design = params$design,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,
        carryover = params$carryover,
        effect_size = NA, se = NA, t_value = NA, p_value = NA, significant = NA
      )
    })

    results <- bind_rows(results, model_result)
  }
}

# ============================================================================
# SUMMARIZE AND SAVE
# ============================================================================

cat("\n", strrep("=", 50), "\n")
cat("RESULTS - EVENLY-SPACED DESIGNS\n")
cat(strrep("=", 50), "\n\n")

summary_results <- results %>%
  group_by(design, biomarker_moderation, biomarker_correlation, carryover) %>%
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(summary_results)

# ============================================================================
# VISUALIZATION
# ============================================================================

library(viridis)

plot_power_heatmap <- function(data) {
  design_labels <- c(
    "ol" = "OL\n(Design 1)",
    "crossover" = "Crossover\n(Design 3)",
    "parallel" = "Parallel"
  )

  data <- data %>%
    mutate(
      carryover = factor(carryover),
      biomarker_moderation = factor(
        biomarker_moderation,
        levels = rev(sort(unique(biomarker_moderation)))
      ),
      design = factor(
        design_labels[as.character(design)],
        levels = design_labels
      )
    )

  p <- ggplot(
    data,
    aes(x = carryover, y = biomarker_moderation, fill = power)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.0f%%", power * 100)),
      color = "black",
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      name = "Power",
      low = "#d73027",
      mid = "#fee08b",
      high = "#1a9850",
      midpoint = 0.5,
      limits = c(0, 1),
      labels = scales::percent
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_wrap(~design, ncol = 3) +
    labs(
      title = "Statistical Power - Evenly-Spaced Designs (4 timepoints)",
      subtitle = sprintf(
        "N=%d participants, %d iterations per condition | Schedule: weeks %s",
        n_participants, n_iterations, paste(measurement_weeks, collapse = ", ")
      ),
      x = "Carryover (proportion retained)",
      y = "Biomarker Moderation\n(effect size)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = NA)
    )

  return(p)
}

p_heatmap <- plot_power_heatmap(summary_results)
print(p_heatmap)

# Save outputs
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "power_heatmap_evenly_spaced.pdf"),
  p_heatmap,
  width = 12,
  height = 7
)

save(
  results, summary_results,
  file = file.path(output_dir, "simulation_evenly_spaced_results.RData")
)

cat("\nDone! Results saved to", output_dir, "\n")
cat("- power_heatmap_evenly_spaced.pdf\n")
cat("- simulation_evenly_spaced_results.RData\n")

# Close log file
sink()
